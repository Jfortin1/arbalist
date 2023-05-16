#include "Rcpp.h"
#include "FragmentParser.hpp"
#include "hdf5_utils.h"

#include <unordered_map>
#include <vector>

struct RegionCounter {
    RegionCounter(Rcpp::CharacterVector seqnames, Rcpp::List region_ids, Rcpp::List region_starts, Rcpp::List region_ends, Rcpp::Nullable<Rcpp::CharacterVector> cellnames) {
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
                Rcpp::IntegerVector(region_ids[i]),
                Rcpp::IntegerVector(region_starts[i]),
                Rcpp::IntegerVector(region_ends[i])
            );
        }

        current_seq_it = seq_to_id.end();
    }

public:
    bool known_cells;
    std::unordered_map<std::string, int> cell_to_id;

    struct Sequence {
        Sequence() = default;
        Sequence(Rcpp::IntegerVector ri, Rcpp::IntegerVector rs, Rcpp::IntegerVector re) : ids(std::move(ri)), starts(std::move(rs)), ends(std::move(re)) {}
        bool visited = false;
        Rcpp::IntegerVector ids, starts, ends;
    };
    std::unordered_map<std::string, Sequence> seq_to_id;

    std::vector<std::vector<int> > collected;
    std::unordered_map<std::string, Sequence>::const_iterator current_seq_it;

public:
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
            auto sIt = seq_to_id.find(seq_name);
            if (sIt == seq_to_id.end()) {
                throw std::runtime_error("unrecognized sequence name '" + seq_name + "' on line " + std::to_string(line_number));
            }
            sIt->second.visited = true;

            current_seq_it = sIt;
            start_region_index = 0;
            end_region_index = 0;

        } else if (current_seq_it->first != seq_name) {
            auto sIt = seq_to_id.find(seq_name);
            if (sIt == seq_to_id.end()) {
                throw std::runtime_error("unrecognized sequence name '" + seq_name + "' on line " + std::to_string(line_number));
            }

            if ((sIt->second).visited) {
                throw std::runtime_error("entries in the fragment file have unordered sequence names (line "  + std::to_string(line_number) + ")");
            }
            sIt->second.visited = true;

            current_seq_it = sIt;
            start_region_index = 0;
            end_region_index = 0;
        }

        if (end_pos <= start_pos) {
            throw std::runtime_error("fragment end (" + std::to_string(end_pos) + ") should be greater than the fragment start (" + std::to_string(start_pos) +") on line " + std::to_string(line_number));
        }

        const auto& ids = (current_seq_it->second).ids;
        const auto& starts = (current_seq_it->second).starts;
        const auto& ends = (current_seq_it->second).ends;

        int nregions = ends.size();
        if (start_region_index == nregions) {
            return;
        }

        // Jumping ahead to the first region that ends after the fragment start.
        // 'start_region_index' should only ever increase, as the fragments are sorted by starts.
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
            // Moving backwards: searching for the region that starts before the fragment end.
            do { 
                --end_region_index;
            } while (end_region_index > start_region_index && starts[end_region_index] > end_pos);

            if (starts[end_region_index] <= end_pos && ends[end_region_index] > end_pos) {
                has_end = true;
                end_id = ids[end_region_index];
            }
        } else {
            // Moving forwards: searching for the region that ends after the fragment end.
            while (end_region_index < nregions && ends[end_region_index] <= end_pos) {
                ++end_region_index;
            }

            if (end_region_index < nregions && starts[end_region_index] <= end_pos) {
                has_end = true;
                end_id = ids[end_region_index];
            }
        }

        // Not counting a fragment twice if it overlaps the same gene.
        bool has_start = starts[start_region_index] <= start_pos;
        int start_id = (has_start ? ids[start_region_index] : -1);

        if (has_start && has_end) {
            collected[cid].push_back(start_id);
            if (start_id != end_id) {
                collected[cid].push_back(end_id);
            }
        } else if (has_start) {
            collected[cid].push_back(start_id);
        } else if (has_end) {
            collected[cid].push_back(end_id);
        }
    }
};

// [[Rcpp::export(rng=false)]]
SEXP fragments_to_regions(
    std::string fragment_file, 
    std::string output_file, 
    std::string output_group, 
    Rcpp::CharacterVector seqnames, 
    Rcpp::List region_ids, 
    Rcpp::List region_starts, 
    Rcpp::List region_ends, 
    Rcpp::Nullable<Rcpp::CharacterVector> cellnames,
    int num_regions,
    int deflate_level,
    int chunk_dim)
{
    RegionCounter counter(seqnames, region_ids, region_starts, region_ends, cellnames);
    {
        FragmentParser parser(&counter);
        parser.run(fragment_file);
    }

    dump_sparse_matrix(output_file, output_group, counter.collected, num_regions, deflate_level, chunk_dim);

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
