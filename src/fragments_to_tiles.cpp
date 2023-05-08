#include "Rcpp.h"
#include "byteme/SomeFileReader.hpp"
#include "H5Cpp.h"

#include <unordered_map>
#include <vector>
#include <queue>

struct TileParser {
    TileParser(int ts, Rcpp::IntegerVector seqlengths, Rcpp::CharacterVector seqnames, Rcpp::Nullable<Rcpp::CharacterVector> cellnames) : tile_size(ts) {
        known_cells = cellnames.isNotNull();
        if (known_cells) {
            Rcpp::CharacterVector my_cellnames(cellnames);
            for (size_t i = 0, end = my_cellnames.size(); i < end; ++i) {
                cell_to_id[Rcpp::as<std::string>(my_cellnames[i])] = i;
            }
            collected.resize(my_cellnames.size());
        }

        if (seqlengths.size() != seqnames.size()) {
            throw std::runtime_error("'seqlengths' and 'seqnames' should be of the same length");
        }

        int offset = 0;
        for (size_t i = 0, end = seqnames.size(); i < end; ++i) {
            seq_to_id[Rcpp::as<std::string>(seqnames[i])] = Sequence(i, seqlengths[i], offset);
            offset += std::ceil(static_cast<double>(seqlengths[i]) / tile_size);
        }

        current_tile_seq_it = seq_to_id.end();
    }

public:
    int tile_size;

    bool known_cells;
    std::unordered_map<std::string, int> cell_to_id;

    struct Sequence {
        Sequence() = default;
        Sequence(int s, int l, int o) : seq_id(s), length(l), offset(o) {}
        int seq_id;
        int length;
        int offset;
    };
    std::unordered_map<std::string, Sequence> seq_to_id;

public:
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
    bool empty = true; 

    std::string seq_name;
    int start_pos = 0;
    int end_pos = 0;
    std::string cell_name;

    std::unordered_map<std::string, Sequence>::const_iterator current_tile_seq_it;
    int last_start_pos = 0;

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
        int cid;
        auto cIt = cell_to_id.find(cell_name);
        if (known_cells) {
            if (cIt == cell_to_id.end()) {
                return; // ignoring unknown cell.
            }
            cid = cIt->second;
        } else {
            if (cIt == cell_to_id.end()) {
                cid = collected.size();
                collected.resize(cid + 1);
                cell_to_id[cell_name] = cid;
            } else {
                cid = cIt->second;
            }
        }

        bool sameseq = true;
        if (current_tile_seq_it == seq_to_id.end()) {
            sameseq = false; // for the start.
            current_tile_seq_it = seq_to_id.find(seq_name);
            if (current_tile_seq_it == seq_to_id.end()) {
                throw std::runtime_error("unrecognized sequence name '" + seq_name + "'");
            }

        } else if (current_tile_seq_it->first != seq_name) {
            sameseq = false;
            auto sIt = seq_to_id.find(seq_name);
            if (sIt == seq_to_id.end()) {
                throw std::runtime_error("unrecognized sequence name '" + seq_name + "'");
            }

            auto sid = (sIt->second).seq_id;
            if ((current_tile_seq_it->second).seq_id > sid) {
                throw std::runtime_error("order of sequences in fragment file differs from that in 'seqlengths'");
            }

            flush_all_leftovers(); // from the previous seqname, if any exists.
            current_tile_seq_it = sIt;
            last_start_pos = start_pos;
        }

        auto slength = (current_tile_seq_it->second).length;
        if (slength <= start_pos || slength <= end_pos) {
            throw std::runtime_error("fragment boundaries out of range of the sequence length");
        }

        int local_tile_id = (start_pos / tile_size);
        auto soffset = (current_tile_seq_it->second).offset;
        int global_tile_id = local_tile_id + soffset;

        if (sameseq) {
            if (start_pos < last_start_pos) {
                throw std::runtime_error("fragment file should be ordered by start position");
            }
            last_start_pos = start_pos;

            // Processing the leftovers up to and including the current tile.
            while (leftovers.size() && leftovers.top().tile_id <= global_tile_id) {
                const auto& top = leftovers.top();
                add_tile(top.cell_id, top.tile_id);
                leftovers.pop();
            }
        }

        add_tile(cid, global_tile_id);

        // Recording the end position.
        if (end_pos <= start_pos) {
            throw std::runtime_error("fragment end should be greater than the fragment start");
        }

        int global_end_tile_id = (end_pos / tile_size) + soffset;
        if (global_end_tile_id == global_tile_id) {
            add_tile(cid, global_tile_id);
        } else {
            leftovers.emplace(cid, global_end_tile_id);
        }
   }

public:
    void parse_buffer(const char* ptr, size_t available) {
        for (size_t i = 0; i < available; ++i) {
            switch(ptr[i]) {
                case '\n':
                    if (field != 4 || empty) {
                        throw std::runtime_error("expected 5 non-empty fields per non-comment line");
                    }
                    add_record();

                    // Setting everything back.
                    field = 0;
                    empty = true;

                    seq_name.clear();
                    start_pos = 0;
                    end_pos = 0;
                    cell_name.clear();
                    break;
                case '\t':
                    if (empty) {
                        throw std::runtime_error("expected non-empty field in the fragment file");
                    }
                    ++field;
                    empty = true;
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
                                throw std::runtime_error("only digits should be present in the end position field");
                            }
                            end_pos *= 10;
                            end_pos += ptr[i] - '0';
                            break;
                        case 3:
                            cell_name += ptr[i];
                            break;
                    }
                    empty = false;
                    break;
            }
        }
    }
};

// [[Rcpp::export(rng=false)]]
Rcpp::List dump_fragments_to_files(
    std::string fragment_file, 
    int tile_size, 
    std::string output_file, 
    std::string output_group, 
    Rcpp::IntegerVector seqlengths, 
    Rcpp::CharacterVector seqnames, 
    Rcpp::Nullable<Rcpp::CharacterVector> cellnames,
    double previous_nonzero) 
{
    TileParser parser(tile_size, seqlengths, seqnames, cellnames);

    {
        byteme::SomeFileReader gzreader(fragment_file);
        bool remaining = false;

        /*** Reading through all the comments at the start. ***/
        bool in_comment = false;
        do {
            remaining = gzreader();
            const unsigned char* buffer = gzreader.buffer();
            size_t available = gzreader.available();
            auto ptr = reinterpret_cast<const char*>(buffer);

            bool breakout = false;
            for (size_t position = 0; position < available; ++position) {
                if (!in_comment) {
                    if (ptr[position] == '#') {
                        in_comment = true;
                    } else {
                        parser.parse_buffer(ptr + position, available - position);
                        breakout = true;
                        break;
                    }
                } else {
                    if (ptr[position] == '\n') {
                        in_comment = false;
                    }
                }
            }

            if (breakout) {
                break;
            }
        } while (remaining);

        /*** Now actually doing the parsing.  ***/
        do {
            remaining = gzreader();
            const unsigned char* buffer = gzreader.buffer();
            size_t available = gzreader.available();

            // Parsing one line at a time and putting it into our store.
            auto ptr = reinterpret_cast<const char*>(buffer);
            parser.parse_buffer(ptr, available);
        } while (remaining);

        if (parser.seq_name.size()) { // add the last record manually, when there's no trailing newline.
            if (parser.field != 4 || parser.empty) {
                throw std::runtime_error("last line does not have 5 fields");
            } else {
                parser.add_record(); 
            }
        }
        parser.flush_all_leftovers();
    }

    // Finally, dumping it all to HDF5.
    H5::H5File fhandle(output_file, H5F_ACC_RDWR);
    H5::Group ghandle = fhandle.openGroup(output_group);

    const auto& collected = parser.collected;
    std::vector<hsize_t> gathered(collected.size());
    hsize_t extension = 0, baseline = previous_nonzero;
    for (size_t i = 0, end = collected.size(); i < end; ++i) {
        size_t n = collected[i].size();
        extension += n;
        gathered[i] = extension + baseline;
    }
    std::cout << "Gathered " << gathered.size() << std::endl;

    {
        H5::DataSet phandle = ghandle.openDataSet("indptr");
        hsize_t sofar = 0;
        phandle.getSpace().getSimpleExtentDims(&sofar);
        hsize_t count = gathered.size();
        hsize_t combined = sofar + gathered.size();
        phandle.extend(&combined);

        H5::DataSpace dataspace(1, &combined), memspace(1, &count);
        dataspace.selectHyperslab(H5S_SELECT_SET, &count, &sofar);
        phandle.write(gathered.data(), H5::PredType::NATIVE_HSIZE, memspace, dataspace);
    }

    H5::DataSet dhandle = ghandle.openDataSet("data");
    hsize_t sofar = 0;
    dhandle.getSpace().getSimpleExtentDims(&sofar);
    hsize_t combined = sofar + extension;
    dhandle.extend(&combined);

    H5::DataSet ihandle = ghandle.openDataSet("indices");
    {
        hsize_t isofar = 0;
        ihandle.getSpace().getSimpleExtentDims(&isofar);
        if (sofar != isofar) {
            throw std::runtime_error("lengths of index and data HDF5 datasets should be the same");
        }
        ihandle.extend(&combined);
    }

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

        dataspace.selectHyperslab(H5S_SELECT_SET, &count, &sofar);
        memspace.setExtentSimple(1, &count);
        dhandle.write(dbuffer.data(), H5::PredType::NATIVE_INT, memspace, dataspace);
        ihandle.write(ibuffer.data(), H5::PredType::NATIVE_INT, memspace, dataspace);

        sofar += count;
    }

    Rcpp::List output(2);
    output[0] = static_cast<double>(extension + baseline);
    if (cellnames.isNotNull()) {
        output[1] = R_NilValue;
    } else {
        Rcpp::CharacterVector names(parser.cell_to_id.size());
        for (const auto& x : parser.cell_to_id) {
            names[x.second] = x.first;
        }
        output[1] = names;
    }
    return output;
}
