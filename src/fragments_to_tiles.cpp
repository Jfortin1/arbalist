#include "Rcpp.h"
#include "byteme/SomeFileReader.hpp"
#include "H5Cpp.h"

#include <unordered_map>
#include <vector>
#include <queue>

/*************************************************
 *** Parsing and counting fragments into tiles ***
 *************************************************/

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

        number_of_bins = 0;
        for (size_t i = 0, end = seqnames.size(); i < end; ++i) {
            seq_to_id[Rcpp::as<std::string>(seqnames[i])] = Sequence(i, seqlengths[i], number_of_bins);
            number_of_bins += std::ceil(static_cast<double>(seqlengths[i]) / tile_size);
        }

        current_tile_seq_it = seq_to_id.end();
    }

public:
    int tile_size;
    int number_of_bins;

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
    bool cycle = false;

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

        if (end_pos <= start_pos) {
            throw std::runtime_error("fragment end should be greater than the fragment start");
        }

        // Avoid auto-correlations between adjacent tiles and over-representation of even counts
        // by only counting one end of each fragment. This avoids the problems from effective
        // correlations when both ends are treated as "independent". Also ensures that the
        // total count for each cell is equal to the number of fragments.
        cycle = !cycle;
        if (cycle) {
            add_tile(cid, global_tile_id);
        } else {
            int global_end_tile_id = (end_pos / tile_size) + soffset;
            if (global_end_tile_id == global_tile_id) {
                add_tile(cid, global_tile_id);
            } else {
                leftovers.emplace(cid, global_end_tile_id);
            }
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

/****************************
 *** Saving HDF5 matrices ***
 ****************************/

inline H5::DataSet create_1d_compressed_hdf5_dataset(H5::Group& location, const H5::DataType& dtype, const std::string& name, hsize_t length, int deflate_level, hsize_t chunk) {
    H5::DataSpace dspace(1, &length);
 	H5::DSetCreatPropList plist;

    if (deflate_level >= 0 && length) {
        plist.setDeflate(deflate_level);
        if (chunk > length) {
            plist.setChunk(1, &length);
        } else {
            plist.setChunk(1, &chunk);
        }
    }

    return location.createDataSet(name, dtype, dspace, plist);
}

/*************************
 *** Exported function ***
 *************************/

// [[Rcpp::export(rng=false)]]
SEXP dump_fragments_to_files(
    std::string fragment_file, 
    int tile_size, 
    std::string output_file, 
    std::string output_group, 
    Rcpp::IntegerVector seqlengths, 
    Rcpp::CharacterVector seqnames, 
    Rcpp::Nullable<Rcpp::CharacterVector> cellnames,
    int deflate_level,
    int chunk_dim)
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
    H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
    fapl.setCache(0, 511, chunk_dim * sizeof(int), 1);
    H5::H5File fhandle(output_file, H5F_ACC_RDWR, H5::FileCreatPropList::DEFAULT, fapl);
    H5::Group ghandle = fhandle.openGroup(output_group);

    const auto& collected = parser.collected;
    std::vector<hsize_t> gathered(collected.size() + 1);
    int max_count = 0;
    for (size_t i = 0, end = collected.size(); i < end; ++i) {
        size_t n = collected[i].size();
        gathered[i + 1] = gathered[i] + n;
        for (const auto& x : collected[i]) {
            if (x.count > max_count) {
                max_count = x.count;
            }
        }
    }

    {
        H5::DataSet phandle = create_1d_compressed_hdf5_dataset(ghandle, H5::PredType::NATIVE_UINT64, "indptr", gathered.size(), deflate_level, chunk_dim);
        hsize_t len = gathered.size();
        H5::DataSpace dataspace(1, &len), memspace(1, &len);
        dataspace.selectAll();
        memspace.selectAll();
        phandle.write(gathered.data(), H5::PredType::NATIVE_HSIZE, memspace, dataspace);
    }

    {
        H5::DataSet shandle = create_1d_compressed_hdf5_dataset(ghandle, H5::PredType::NATIVE_UINT64, "shape", 2, deflate_level, chunk_dim);
        hsize_t len = 2;
        H5::DataSpace dataspace(1, &len), memspace(1, &len);
        dataspace.selectAll();
        memspace.selectAll();
        hsize_t shape[2];
        shape[0] = parser.number_of_bins;
        shape[1] = collected.size();
        shandle.write(shape, H5::PredType::NATIVE_HSIZE, memspace, dataspace);
    }

    H5::DataSet ihandle = create_1d_compressed_hdf5_dataset(ghandle, H5::PredType::NATIVE_UINT32, "indices", gathered.back(), deflate_level, chunk_dim);

    const H5::PredType* dtype;
    if (max_count <= 255) {
        dtype = &(H5::PredType::NATIVE_UINT8);
    } else if (max_count <= 65535) {
        dtype = &(H5::PredType::NATIVE_UINT16);
    } else {
        dtype = &(H5::PredType::NATIVE_UINT32);
    }
    H5::DataSet dhandle = create_1d_compressed_hdf5_dataset(ghandle, *dtype, "data", gathered.back(), deflate_level, chunk_dim);

    std::vector<int> ibuffer;
    std::vector<int> dbuffer;
    H5::DataSpace dataspace(1, &(gathered.back())), memspace;
    hsize_t sofar = 0;

    for (const auto& x : collected) {
        ibuffer.clear();
        dbuffer.clear();

        hsize_t count = x.size();
        ibuffer.reserve(count);
        dbuffer.reserve(count);

        for (const auto& y : x) {
            ibuffer.push_back(y.tile_id);
            dbuffer.push_back(y.count); // (y.count + 1)/2); // Avoid bias towards even counts by dividing by 2, see https://www.biorxiv.org/content/10.1101/2022.05.04.490536v1.full.
        }

        dataspace.selectHyperslab(H5S_SELECT_SET, &count, &sofar);
        memspace.setExtentSimple(1, &count);
        memspace.selectAll();
        dhandle.write(dbuffer.data(), H5::PredType::NATIVE_INT, memspace, dataspace);
        ihandle.write(ibuffer.data(), H5::PredType::NATIVE_INT, memspace, dataspace);

        sofar += count;
    }

    if (cellnames.isNotNull()) {
        return R_NilValue;
    } else {
        Rcpp::CharacterVector names(parser.cell_to_id.size());
        for (const auto& x : parser.cell_to_id) {
            names[x.second] = x.first;
        }
        return names;
    }
}
