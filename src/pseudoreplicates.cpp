#include "Rcpp.h"
#include "FragmentParser.hpp"
#include <iostream>
#include <fstream>

#include <unordered_map>
#include <vector>

struct PseudoBulkMaker {
  PseudoBulkMaker(std::string output_file, Rcpp::Nullable<Rcpp::CharacterVector> cellnames) {
   known_cells = cellnames.isNotNull();
    if (known_cells) {
      Rcpp::CharacterVector my_cellnames(cellnames);
      for (size_t i = 0, end = my_cellnames.size(); i < end; ++i) {
        cell_to_id[Rcpp::as<std::string>(my_cellnames[i])] = i;
      }
      collected.resize(my_cellnames.size());
    }
    output_pseudobulk_file = output_file;
    outfile.open(output_file);
  }
  
public:
  bool known_cells;
  std::unordered_map<std::string, int> cell_to_id;
  std::vector<std::vector<int> > collected;
  std::string output_pseudobulk_file;
  std::ofstream outfile;

public:
  void add(const std::string& seq_name, int start_pos, int end_pos, const std::string& cell_name, int, size_t line_number) {
      auto cIt = cell_to_id.find(cell_name);
      if (cIt != cell_to_id.end()) {
        outfile << seq_name << "\t" << start_pos << "\t" << end_pos << "\t" << cell_name << "\n";
      }
  }
};

// [[Rcpp::export(rng=false)]]
void create_pseudobulk_file(
    Rcpp::Nullable<Rcpp::CharacterVector> fragment_files, 
    std::string output_file, 
    Rcpp::Nullable<Rcpp::CharacterVector> cellnames)
{

  Rcpp::CharacterVector my_fragment_files(fragment_files);
  
  PseudoBulkMaker counter(output_file, cellnames);
  {
    FragmentParser parser(&counter);
    for (size_t i = 0, end = my_fragment_files.size(); i < end; ++i) {
      parser.run(Rcpp::as<std::string>(my_fragment_files[i]));
    }
  }
  
  counter.outfile.close();
  
}