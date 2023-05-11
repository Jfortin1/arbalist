#include "Rcpp.h"
#include "FragmentParser.hpp"
#include "hdf5_utils.hpp"

#include <unordered_map>
#include <vector>
#include <queue>

struct RegionCounter {
    RegionCounter(int ts, Rcpp::CharacterVector seqnames, Rcpp::List region_ids, Rcpp::List region_starts, Rcpp::List region_ends, Rcpp::Nullable<Rcpp::CharacterVector> cellnames) : tile_size(ts) {
        known_cells = cellnames.isNotNull();
        if (known_cells) {
            Rcpp::CharacterVector my_cellnames(cellnames);
            for (size_t i = 0, end = my_cellnames.size(); i < end; ++i) {
                cell_to_id[Rcpp::as<std::string>(my_cellnames[i])] = i;
            }
            collected.resize(my_cellnames.size());
        }

        for (size_t i = 0, end = seqnames.size(); i < end; ++i) {
            seq_to_id[Rcpp::as<std::string>(seqnames[i])] = Sequence(
                i,
                Rcpp::IntegerVector(region_ids[i]),
                Rcpp::IntegerVector(region_starts[i]),
                Rcpp::IntegerVector(region_ends[i])
            );
        }

        current_tile_seq_it = seq_to_id.end();
    }

public:
    bool known_cells;
    std::unordered_map<std::string, int> cell_to_id;

    struct Sequence {
        Sequence() = default;
        Sequence(int id, Rcpp::IntegerVector ri, Rcpp::IntegerVector rs, Rcpp::IntegerVector re) : 
            seq_id(id), ids(std::move(ri)), starts(std::move(rs)), ends(std::move(re)) {}

        int seq_id;
        Rcpp::IntegerVector ids, starts, ends;
    };
    std::unordered_map<std::string, Sequence> seq_to_id;

    std::vector<std::vector<int> > collected;
    bool cycle = false;
    std::unordered_map<std::string, Sequence>::const_iterator current_seq_it;

public:
    struct EndPosition {
        EndPosition(int c, int s, int e) : cell_id(c), start_id(s), end_pos(e) {}
        int cell_id;
        int start_id;
        int end_pos;
    };

    struct CompareEndPosition {
        bool operator()(const EndPosition& left, const EndPosition& right) const {
            return left.end_pos > right.end_pos;
        }
    }

    std::priority_queue<Ends, std::vector<Ends>, CompareEndPosition> end_positions;
    int start_region_index = 0, end_region_index = 0;

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

        if (current_seq_it == seq_to_id.end()) {
            current_seq_it = seq_to_id.find(seq_name);
            if (current_seq_it == seq_to_id.end()) {
                throw std::runtime_error("unrecognized sequence name '" + seq_name + "' on line " + std::to_string(line_number));
            }
            start_region_index = 0;
            end_region_index = 0;

        } else if (current_seq_it->first != seq_name) {
            auto sIt = seq_to_id.find(seq_name);
            if (sIt == seq_to_id.end()) {
                throw std::runtime_error("unrecognized sequence name '" + seq_name + "' on line " + std::to_string(line_number));
            }

            auto sid = (sIt->second).seq_id;
            if ((current_seq_it->second).seq_id > sid) {
                throw std::runtime_error("order of sequences in fragment file differs from that in 'seqlengths' (line "  + std::to_string(line_number) + ")");
            }

            current_seq_it = sIt;
            start_region_index = 0;
            end_region_index = 0;
        }

        if (end_pos <= start_pos) {
            throw std::runtime_error("fragment end (" + std::to_string(end_pos) + ") should be greater than the fragment start (" + std::to_string(start_pos) +") on line " + std::to_string(line_number));
        }

        const auto& ids = current_seq_it->ids;
        const auto& starts = current_seq_it->starts;
        const auto& ends = current_seq_it->ends;

        int nregions = ends.size();
        if (start_region_index == nregions) {
            return;
        }

        if (ends[start_region_index] <= start_pos) {
            do {
                ++start_region_index;
            } while (start_region_index < nregions && ends[start_region_index] <= start_pos);

            end_region_index = start_region_index;
            if (start_region_index == nregions) {
                return;
            }
        }

        bool has_end = false;
        int end_id = -1;

        if (end_region_index == nregions || starts[end_region_index] > end_pos) {
            // Walking backwards.
            do { 
                --end_region_index;
            } while (end_region_index > start_region_index && starts[end_region_index] > end_pos);

            if (ends[end_region_index] > end_pos) {
                has_end = true;
                end_id = ends[end_region_index];
            }
        } else {
            // Walking forwards.
            while (end_region_index < nregions && ends[end_region_index] <= end_pos) {
                ++end_region_index;
            }

            if (end_region_index < nregions && starts[end_region_index] <= end_pos) {
                has_end = true;
                end_id = ends[end_region_index];
            }
        }

        bool has_start = starts[start_region_index] <= start_pos;
        int start_id = (has_? ids[start_region_index] : -1);

        if (has_start && has_end) {
            if (start_id == end_id) {
                collected[cid].push_back(start_id);
            }
        } else if (has_start) {
            collected[cid].push_back(start_id);
        } else if (has_end) {
            collected[cid].push_back(end_id);
        }
    }
};

// [[Rcpp::export(rng=false)]]
SEXP dump_fragments_to_files(
    std::string fragment_file, 
    int tile_size, 
    std::string output_file, 
    std::string output_group, 
    Rcpp::CharacterVector seqnames, 
    Rcpp::List region_ids, 
    Rcpp::List region_starts, 
    Rcpp::List region_ends, 
    Rcpp::Nullable<Rcpp::CharacterVector> cellnames,
    int deflate_level,
    int chunk_dim)
{
    RegionCounter counter(tile_size, seqlengths, seqnames, cellnames);
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
