#include "Rcpp.h"
#include "H5Cpp.h"
#include "FragmentParser.hpp"

#include <unordered_map>
#include <vector>
#include <queue>

struct TileCounter {
    TileCounter(int ts, Rcpp::IntegerVector seqlengths, Rcpp::CharacterVector seqnames, Rcpp::Nullable<Rcpp::CharacterVector> cellnames) : tile_size(ts) {
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

    std::vector<std::vector<int> > collected;
    bool cycle = false;
    std::unordered_map<std::string, Sequence>::const_iterator current_tile_seq_it;

public:
    void add(const std::string& seq_name, int start_pos, int end_pos, const std::string& cell_name, int, size_t line_number) {
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

        if (current_tile_seq_it == seq_to_id.end()) {
            current_tile_seq_it = seq_to_id.find(seq_name);
            if (current_tile_seq_it == seq_to_id.end()) {
                throw std::runtime_error("unrecognized sequence name '" + seq_name + "' on line " + std::to_string(line_number));
            }

        } else if (current_tile_seq_it->first != seq_name) {
            auto sIt = seq_to_id.find(seq_name);
            if (sIt == seq_to_id.end()) {
                throw std::runtime_error("unrecognized sequence name '" + seq_name + "' on line " + std::to_string(line_number));
            }

            auto sid = (sIt->second).seq_id;
            if ((current_tile_seq_it->second).seq_id > sid) {
                throw std::runtime_error("order of sequences in fragment file differs from that in 'seqlengths' (line "  + std::to_string(line_number) + ")");
            }

            current_tile_seq_it = sIt;
        }

        auto slength = (current_tile_seq_it->second).length;
        if (slength <= start_pos || slength <= end_pos) {
            throw std::runtime_error("fragment boundaries (" + std::to_string(start_pos) + ":" + std::to_string(end_pos) + ") out of range of the sequence length on line " + std::to_string(line_number));
        }

        if (end_pos <= start_pos) {
            throw std::runtime_error("fragment end (" + std::to_string(end_pos) + ") should be greater than the fragment start (" + std::to_string(start_pos) +") on line " + std::to_string(line_number));
        }

        // Avoid auto-correlations between adjacent tiles and over-representation of even counts
        // by only counting one end of each fragment. This avoids the problems from effective
        // correlations when both ends are treated as "independent". Also ensures that the
        // total count for each cell is equal to the number of fragments.
        cycle = !cycle;
        int global_tile_id = (current_tile_seq_it->second).offset;
        if (cycle) {
            global_tile_id += (start_pos / tile_size);
        } else {
            global_tile_id += (end_pos / tile_size);
        }
        collected[cid].push_back(global_tile_id);
    }
};

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
    TileCounter counter(tile_size, seqlengths, seqnames, cellnames);
    {
        FragmentParser parser(&counter);
        parser.run(fragment_file);
    }

    // Dumping it all to HDF5.
    H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
    fapl.setCache(0, 511, chunk_dim * 10 * sizeof(int), 1);
    H5::H5File fhandle(output_file, H5F_ACC_RDWR, H5::FileCreatPropList::DEFAULT, fapl);
    H5::Group ghandle = fhandle.openGroup(output_group);

    auto& collected = counter.collected;
    std::vector<hsize_t> gathered(collected.size() + 1);
    int max_count = 0;
    for (size_t i = 0, end = collected.size(); i < end; ++i) {
        auto& x = collected[i];
        std::sort(x.begin(), x.end());

        int nunique = 0;
        int last = -1;
        int count = 0;
        for (auto y : x) {
            if (y > last) {
                ++nunique;
                if (count > max_count) {
                    max_count = count;
                }
                count = 1;
                last = y;
            } else {
                ++count;
            }
        }

        gathered[i + 1] = gathered[i] + nunique;
        if (count > max_count) {
            max_count = count;
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
        shape[0] = counter.number_of_bins;
        shape[1] = collected.size();
        shandle.write(shape, H5::PredType::NATIVE_HSIZE, memspace, dataspace);
    }

    hsize_t total = gathered.back();
    H5::DataSet ihandle = create_1d_compressed_hdf5_dataset(ghandle, H5::PredType::NATIVE_UINT32, "indices", total, deflate_level, chunk_dim);

    const H5::PredType* dtype;
    if (max_count <= 255) {
        dtype = &(H5::PredType::NATIVE_UINT8);
    } else if (max_count <= 65535) {
        dtype = &(H5::PredType::NATIVE_UINT16);
    } else {
        dtype = &(H5::PredType::NATIVE_UINT32);
    }
    H5::DataSet dhandle = create_1d_compressed_hdf5_dataset(ghandle, *dtype, "data", total, deflate_level, chunk_dim);

    std::vector<int> ibuffer;
    std::vector<int> dbuffer;
    H5::DataSpace dataspace(1, &total), memspace;
    hsize_t sofar = 0;

    for (auto& x : collected) {
        ibuffer.clear();
        dbuffer.clear();

        // Saving the damn thing in reverse order, like an idiot.
        for (auto y : x) {
            if (ibuffer.empty() || y != ibuffer.back()) {
                ibuffer.push_back(y);
                dbuffer.push_back(1);
            } else {
                ++(dbuffer.back());
            }
        }

        hsize_t count = ibuffer.size();
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
        Rcpp::CharacterVector names(counter.cell_to_id.size());
        for (const auto& x : counter.cell_to_id) {
            names[x.second] = x.first;
        }
        return names;
    }
}
