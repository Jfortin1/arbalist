#include "Rcpp.h"
#include "FragmentParser.hpp"
#include "hdf5_utils.h"

#include <vector>

struct SizeCounter {
  SizeCounter(int user_max_size) {
    max_size_counted = user_max_size;
    fragment_size_counts.reserve(max_size_counted);
    for (size_t i = 0, end = max_size_counted; i < end; ++i) {
      fragment_size_counts[i] = 0;
    }
  }
  
public:
  std::vector<int> fragment_size_counts;
  int max_size_counted;
  int max_size_seen = 0;
  
public:
  void add(const std::string& seq_name, int start_pos, int end_pos, const std::string& cell_name, int, size_t line_number) {
    int fragment_size = end_pos - start_pos-1;
    if(fragment_size < max_size_counted) {
      fragment_size_counts[fragment_size] = fragment_size_counts[fragment_size] + 1;
    } else {
      if(fragment_size > max_size_seen) {
        max_size_seen = fragment_size;
      }
    }
    return;
  }
};

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector count_fragment_size_distributions(std::string fragment_file) 
{
  SizeCounter counter(5500);
  {
    FragmentParser parser(&counter);
    parser.run(fragment_file);
  }
  
  if(counter.max_size_counted < counter.max_size_seen) {
    std::printf("Warning: Observed fragment sizes larger than the distribution we are recording. Max fragment size observed is %i\n", counter.max_size_seen);
  }
  
  Rcpp::IntegerVector res(counter.max_size_counted);
  
  for (size_t i = 0, end = counter.max_size_counted; i < end; ++i) {
    res[i] = counter.fragment_size_counts[i];
  }
  
  return res;
  
}
