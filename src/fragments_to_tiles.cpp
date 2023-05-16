#include "Rcpp.h"
#include "H5Cpp.h"
#include "FragmentParser.hpp"
#include "hdf5_utils.h"

#include <unordered_map>
#include <vector>

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

        int offset = (current_tile_seq_it->second).offset;
        int start_id = (start_pos / tile_size) + offset;
        collected[cid].push_back(start_id);

        int end_id = (end_pos / tile_size) + offset;
        if (start_id != end_id) {
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

    dump_sparse_matrix(output_file, output_group, counter.collected, counter.number_of_bins, deflate_level, chunk_dim);

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
