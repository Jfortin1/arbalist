#ifndef HDF5_UTILS_HPP
#define HDF5_UTILS_HPP

#include <vector>
#include <string>

void dump_sparse_matrix(const std::string&, const std::string&, std::vector<std::vector<int> >&, int, int, int);

#endif
