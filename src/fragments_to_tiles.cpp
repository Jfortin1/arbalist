#include "Rcpp.h"
#include "byteme/GzipFileReader.hpp"

#include <unordered_map>
#include <vector>
#include <queue>

class TileParser {
    TileParser(int ts, Rcpp::IntegerVector seqlengths, Rcpp::CharacterVector seqnames, Rcpp::CharacterVector cellnames) : tile_size(ts) {
        for (size_t i = 0, end = cellnames.size(); i < end; ++i) {
            cell_to_id[Rcpp::as<std::string>(cellnames[i])] = i;
        }
        collected.resize(cellnames.size());

        if (seqlengths.size() != seqnames.size()) {
            throw std::runtime_error("'seqlengths' and 'seqnames' should be of the same length");
        }
        for (size_t i = 0, end = cellnames.size(); i < end; ++i) {
            seq_to_id[Rcpp::as<std::string>(seqnames[i])] = std::pair<int, int>(i, std::ceil(static_cast<double>(seqlengths[i]) / tile_size));
        }
    }

public:
    int tile_size;
    std::unordered_map<std::string, int> cell_to_id;
    std::unordered_map<std::string, std::pair<int, int> > seq_to_id;

    struct Position {
        Position(int c, int t) : cell_id(c), tile_id(t) {}
        int cell_id;
        int tile_id;
    };

    struct Compare {
        // Returns true if left > right, so that the priority queue
        // sorts in increasing order.
        bool operator()(const Position& left, const Position& right) const {
            return left.tile_id > right.tile_id;
        }
    };

    std::priority_queue<Position, std::vector<Position>, Compare> leftovers;

    struct Tile {
        Tile(int t, int c) : tile_id(t), count(c) {}
        int tile_id;
        int count;
    };
    std::vector<std::vector<Tile> > collected;

public:
    unsigned char field = 0;
    std::string seq_name;
    int start_pos = 0;
    int end_pos = 0;
    std::string cell_name;

    void add_tile(int cell_id, int tile_id) {
        auto& vec = collected[cell_id];
        if (vec.empty() || vec.back().tile_id != tile_id) {
            vec.emplace_back(tile_id, 1);
        } else {
            ++(vec.back().count);
        }
    }

    void flush_all_leftovers() {
        while (leftovers.size()) {
            const auto& top = leftovers.top();
            add_tile(top.cell_id, top.tile_id);
            leftovers.pop();
        }
    }

    void add_record() {
        auto cIt = cellnames.find(cell_name);
        if (cIt == cellnames.end()) {
            throw std::runtime_error("unrecognized cell name '" + cell_name + "' in fragment file '" + fragment_file + "'");
        }

        auto sIt = seqnames.find(seq_name);
        if (sIt == seqnames.end()) {
            throw std::runtime_error("unrecognized sequence name '" + seq_name + "' in fragment file '" + fragment_file + "'");
        }
        auto sid = (sIt->second).first;
        auto soffset = (sIt->second).second;

        int local_tile_id = (start_pos / tile_size);
        int global_tile_id = local_tile_id + soffset;

        // Adding the start position.
        if (current_tile_seq_id != sid) {
            flush_all_leftovers(); // from the previous seqname, if any exists.
            current_tile_seq_id = sid;
            current_tile_start = local_tile_id * tile_size;
            current_tile_end = current_tile_start + tile_size;

        } else {
            if (start_pos < current_tile_start) {
                throw std::runtime_error("fragment file should be ordered by start position");
            } 

            if (start_pos < current_tile_end) {
                add_tile(cIt->second, global_tile_id);
            } else {
                // Processing the leftovers up to and including the current tile.
                while (leftovers.size() && leftovers.top().tile_id <= global_tile_id) {
                    const auto& top = leftovers.top();
                    add_tile(top.cell_id, top.tile_id);
                    leftovers.pop();
                }

                add_tile(cIt->second, global_tile_id);
                current_tile_start = local_tile_id * tile_size;
                current_tile_end = current_tile_start + tile_size;
            }
        }

        // Recording the end position.
        if (end_pos <= start_pos) {
            throw std::runtime_error("fragment end should be greater than the fragment start");
        }
        if (end_pos < current_tile_end) {
            add_tile(cIt->second, global_tile_id);
        } else {
            leftovers.emplace_back(cIt->second, (end_pos / tile_size) + soffset);
        }
   }

public:
    void parse_buffer(const char* ptr, size_t available) {
        for (size_t i = 0; i < available; ++i) {
            switch(ptr[i]) {
                case '\n':
                    if (field != 4) {
                        throw std::runtime_error("expected 4 fields per non-comment line");
                    }
                    add_record();

                    // Setting everything back.
                    field = 0;
                    seq_name.clear();
                    start_pos = 0;
                    end_pos = 0;
                    cell_name.clear();
                    break;
                case '\t':
                    ++field;
                    break;
                default:
                    switch (field) {
                        case 0:
                            seq_name += ptr[i];
                            break;
                        case 1:
                            if (!std::isdigit(ptr[i])) {
                                throw std::runtime_error("only digits should be present in the start position field");
                            }
                            start_pos *= 10;
                            start_pos += ptr[i] - '0';
                            break;
                        case 2:
                            if (!std::isdigit(ptr[i])) {
                                throw std::runtime_error("only digits should be present in the start position field");
                            }
                            end_pos *= 10;
                            end_pos += ptr[i] - '0';
                            break;
                        case 3:
                            cell_name += ptr[i];
                            break;
                    }
                    break;
            }
        }
    }
};

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector dump_fragments_to_files(
    std::string fragment_file, 
    int tile_size, 
    std::string output_file, 
    std::string output_group, 
    Rcpp::NumericVector seqlengths, 
    Rcpp::CharacterVector seqnames, 
    Rcpp::CharacterVector cellnames) 
{
    TileParser parser(tile_size, seqlengths, seqnames, cellnames);

    {
        byteme::GzipFileReader gzreader(fragment_file);
        bool remaining = false;

        /*** Reading through all the comments at the start. ***/
        bool in_comment = false;
        do {
            remaining = gzreader();
            const unsigned char* buffer = gzreader.buffer();
            size_t available = gzreader.available();
            auto ptr = reinterpret_cast<const char*>(buffer);

            for (size_t position = 0; position < available; ++position) {
                if (!in_comment) {
                    if (ptr[position] == '#') {
                        in_comment = true;
                    } else {
                        parser.parse_buffer(ptr + position, available - position);
                        break;
                    }
                } else {
                    if (ptr[position] == '\n') {
                        in_comment = false;
                    }
                }
            }
        } while (remaining && in_comment);

        /*** Now actually doing the parsing.  ***/
        do {
            remaining = gzreader();
            const unsigned char* buffer = gzreader.buffer();
            size_t available = gzreader.available();

            // Parsing one line at a time and putting it into our store.
            auto ptr = reinterpret_cast<const char*>(buffer);
            parser.parse_buffer(ptr, available);
        } while (remaining);

        if (seqname) {
            parser.add_record(); // add the last record if there's no trailing newline.
        }
        parser.flush_leftovers();
    }

    // Finally, dumping it all to HDF5.
    H5::H5File fhandle(output_file);
    H5::Group ghandle = fhandle.openGroup(output_name);
    H5::DataSet dhandle = ghandle.openDataSet("data");
    H5::DataSet ihandle = ghandle.openDataSet("indices");

    const auto& collected = parser.collected;
    Rcpp::IntegerVector gathered(collected.size());
    size_t extension = 0;
    for (size_t i = 0, end = collected.size(); i < end; ++i) {
        gathered[i] = x.size();
        extension += x.size();
    }

    hsize_t sofar = 0;
    dhandle.getSpace().getSimpleExtentDims(&sofar);
    hsize_t combined = sofar + extension;
    dhandle.extend(combined);

    hsize_t isofar = 0;
    ihandle.getSpace().getSimpleExtentDims(&isofar);
    if (sofar != isofar) {
        throw std::runtime_error("lengths of index and data HDF5 datasets should be the same");
    }
    ihandle.extend(combined);

    std::vector<int> ibuffer;
    std::vector<int> dbuffer;
    H5::DataSpace dataspace(1, &combined), memspace;

    for (const auto& x : collected) {
        ibuffer.clear();
        dbuffer.clear();

        hsize_t count = x.size();
        ibuffer.reserve(count);
        dbuffer.reserve(count);

        for (const auto& y : x) {
            ibuffer.push_back(y.tile_id);
            dbuffer.push_back(y.count);
        }

        dataspace.setHyperslab(H5S_SELECT_SET, &count, &sofar);
        memspace.setExtentSimple(1, &count);
        dhandle.write(dbuffer.data(), H5::PredType::NATIVE_INT, memspace, dataspace);
        ihandle.write(ibuffer.data(), H5::PredType::NATIVE_INT, memspace, dataspace);
    }

    return gathered;
}
