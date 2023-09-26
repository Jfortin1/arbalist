#include "config.h"

#include "Rcpp.h"
#include "irlba/irlba.hpp"
#include "irlba/parallel.hpp"

#include <vector>
#include <algorithm>
#include <cmath>
#include <thread>

// TODO: move back to libscran.
struct TatamiWrapper {
  TatamiWrapper(const tatami::NumericMatrix* m, int n) : mat(m), nthreads(n) {}
  
private:
  const tatami::NumericMatrix* mat;
  int nthreads;
  
public:
  size_t rows() const {
    return mat->nrow();
  }
  
  size_t cols() const {
    return mat->ncol();
  }
  
public:
  struct Workspace {
  private:
    template<bool row_, bool sparse_>
    void initialize(
        const tatami::NumericMatrix* mat, 
        int nthreads, 
        bool direct, 
        std::vector<std::unique_ptr<tatami::Extractor<tatami::DimensionSelectionType::FULL, sparse_, double, int> > >& work,
        std::vector<std::unique_ptr<tatami::Extractor<tatami::DimensionSelectionType::BLOCK, sparse_, double, int> > >& blockwork)
    {
      tatami::Options opt;
      opt.cache_for_reuse = true;
      int mydim = (row_ ? mat->nrow() : mat->ncol());
      int otherdim = (row_ ? mat->ncol() : mat->nrow());
      
      values.reserve(nthreads);
      if constexpr(sparse_) {
        indices.reserve(nthreads);
      }
      
      if (direct || nthreads == 1) {
        int start = 0;
        int per_thread = std::ceil(static_cast<double>(mydim) / nthreads);
        work.reserve(nthreads);
        ranges.reserve(nthreads);
        
        for (int t = 0; t < nthreads; ++t) {
          int start = per_thread * t;
          int end = std::min(mydim, start + per_thread);
          if (end <= start) {
            break;
          }
          
          work.push_back(tatami::new_extractor<row_, sparse_>(mat, opt));
          values.emplace_back(otherdim);
          if constexpr(sparse_) {
            indices.emplace_back(otherdim);
          }
          
          ranges.emplace_back(start, end);
        }
        
      } else {
        int start = 0;
        int per_thread = std::ceil(static_cast<double>(otherdim) / nthreads);
        blockwork.reserve(nthreads);
        staging.reserve(nthreads);
        
        for (int t = 0; t < nthreads; ++t) {
          int end = std::min(otherdim, start + per_thread);
          if (end <= start) {
            break;
          }
          
          int length = end - start;
          blockwork.push_back(tatami::new_extractor<row_, sparse_>(mat, start, length, opt));
          
          values.emplace_back(length);
          if constexpr(sparse_) {
            indices.emplace_back(length);
          }
          
          staging.emplace_back(otherdim);
          start = end;
        }
      }
      
      use_serial = (values.size() == 1);
    }
    
  public:
    Workspace(const tatami::NumericMatrix* mat, int nthreads, bool multiply) { 
      if (mat->prefer_rows()) {
        if (mat->sparse()) {
          initialize<true, true>(mat, nthreads, multiply, swork, sbwork);
        } else {
          initialize<true, false>(mat, nthreads, multiply, dwork, dbwork);
        }
      } else {
        if (mat->sparse()) {
          initialize<false, true>(mat, nthreads, !multiply, swork, sbwork);
        } else {
          initialize<false, false>(mat, nthreads, !multiply, dwork, dbwork);
        }
      }
    }
    
  private:
    bool use_serial;
    
    std::vector<std::vector<double> > values;
    std::vector<std::vector<int> > indices;
    
    std::vector<std::unique_ptr<tatami::Extractor<tatami::DimensionSelectionType::FULL, true, double, int> > > swork;
    std::vector<std::unique_ptr<tatami::Extractor<tatami::DimensionSelectionType::FULL, false, double, int> > > dwork;
    std::vector<std::pair<int, int> > ranges;
    
    std::vector<std::unique_ptr<tatami::Extractor<tatami::DimensionSelectionType::BLOCK, true, double, int> > > sbwork;
    std::vector<std::unique_ptr<tatami::Extractor<tatami::DimensionSelectionType::BLOCK, false, double, int> > > dbwork;
    std::vector<std::vector<double> > staging; // avoid false sharing during parallelized addition.
    
  private:
    template<bool sparse_>
    auto& get_work() {
      if constexpr(sparse_) {
        return swork;
      } else {
        return dwork;
      }
    }
    
    template<bool sparse_>
    auto& get_block_work() {
      if constexpr(sparse_) {
        return sbwork;
      } else {
        return dbwork;
      }
    }
    
  private:
    template<bool sparse_, class Right_>
    void running_serial(int mydim, const Right_& rhs, Eigen::VectorXd& output) {
      auto& current = get_work<sparse_>().front();
      auto vbuffer = values.front().data();
      auto ibuffer = (sparse_ ? indices.front().data() : static_cast<int*>(NULL));
      
      output.setZero();
      int other = current->full_length;
      
      for (int x = 0; x < mydim; ++x) {
        auto mult = rhs.coeff(x);
        
        if constexpr(sparse_) {
          auto found = current->fetch(x, vbuffer, ibuffer);
          for (int i = 0; i < found.number; ++i) {
            output.coeffRef(found.index[i]) += found.value[i] * mult;
          }
        } else {
          auto found = current->fetch(x, vbuffer);
          for (int i = 0; i < other; ++i) {
            output.coeffRef(i) += found[i] * mult;
          }
        }
      }
    }
    
    template<bool sparse_, class Right_>
    void running_parallel(int mydim, const Right_& rhs, Eigen::VectorXd& output) {
      auto& currentwork = get_block_work<sparse_>();
      output.setZero();
      
      IRLBA_CUSTOM_PARALLEL(currentwork.size(), [&](int thread) -> void {
        auto& current = currentwork[thread];
        auto vbuffer = values[thread].data();
        auto ibuffer = (sparse_ ? indices[thread].data() : static_cast<int*>(NULL));
        
        auto& sums = staging[thread]; // do the sums in a staging area to avoid false sharing problems.
        std::fill(sums.begin(), sums.end(), 0);
        auto blockstart = current->block_start;
        auto blocklen = current->block_length;
        
        for (int x = 0; x < mydim; ++x) {
          auto mult = rhs.coeff(x);
          
          if constexpr(sparse_) {
            auto found = current->fetch(x, vbuffer, ibuffer);
            for (int i = 0; i < found.number; ++i) {
              sums[found.index[i]] += found.value[i] * mult;
            }
          } else {
            auto found = current->fetch(x, vbuffer);
            for (int i = 0; i < blocklen; ++i) {
              sums[blockstart + i] += found[i] * mult;
            }
          }
        }
        
        std::copy_n(sums.begin() + blockstart, blocklen, output.begin() + blockstart);
      });
    }
    
  public:
    template<bool sparse_, class Right_>
    void running(int mydim, const Right_& rhs, Eigen::VectorXd& output) {
      if (use_serial) {
        running_serial<sparse_>(mydim, rhs, output);
      } else {
        running_parallel<sparse_>(mydim, rhs, output);
      }
    }
    
  private:
    template<bool sparse_, class Right_, class Extractor_>
    double direct_raw(const Right_& rhs, int x, double* vbuffer, int* ibuffer, Extractor_& current) {
      double total = 0;
      
      if constexpr(sparse_) {
        auto found = current->fetch(x, vbuffer, ibuffer);
        for (int i = 0; i < found.number; ++i) {
          total += found.value[i] * rhs.coeff(found.index[i]);
        }
      } else {
        auto found = current->fetch(x, vbuffer);
        int otherdim = current->full_length;
        for (int i = 0; i < otherdim; ++i) {
          total += found[i] * rhs.coeff(i);
        }
      }
      
      return total;
    }
    
    template<bool sparse_, class Right_>
    void direct_serial(int mydim, const Right_& rhs, Eigen::VectorXd& output) {
      auto& current = get_work<sparse_>().front();
      auto vbuffer = values.front().data();
      auto ibuffer = (sparse_ ? indices.front().data() : static_cast<int*>(NULL));
      
      output.setZero();
      for (int x = 0; x < mydim; ++x) {
        output.coeffRef(x) = direct_raw<sparse_>(rhs, x, vbuffer, ibuffer, current);
      }
    }
    
    template<bool sparse_, class Right_>
    void direct_parallel(int mydim, const Right_& rhs, Eigen::VectorXd& output) {
      auto& currentwork = get_work<sparse_>();
      output.setZero();
      
      IRLBA_CUSTOM_PARALLEL(currentwork.size(), [&](int thread) -> void {
        auto range = ranges[thread];
        auto start = range.first;
        auto end = range.second;
        
        auto& current = currentwork[thread];
        auto vbuffer = values[thread].data();
        auto ibuffer = (sparse_ ? indices[thread].data() : static_cast<int*>(NULL));
        
        for (int x = start; x < end; ++x) {
          output.coeffRef(x) = direct_raw<sparse_>(rhs, x, vbuffer, ibuffer, current);
        }
      });
    }
    
  public:
    template<bool sparse_, class Right_>
    void direct(int mydim, const Right_& rhs, Eigen::VectorXd& output) {
      if (use_serial) {
        direct_serial<sparse_>(mydim, rhs, output);
      } else {
        direct_parallel<sparse_>(mydim, rhs, output); 
      }
    }
  };
  
public:
  Workspace workspace() const {
    return Workspace(mat, nthreads, true);
  }
  
  Workspace adjoint_workspace() const {
    return Workspace(mat, nthreads, false);
  }
  
  template<class Right>
  void multiply(const Right& rhs, Workspace& work, Eigen::VectorXd& output) const {
    if (mat->prefer_rows()) {
      if (mat->sparse()) {
        work.template direct<true>(mat->nrow(), rhs, output);
      } else {
        work.template direct<false>(mat->nrow(), rhs, output);
      }
    } else {
      if (mat->sparse()) {
        work.template running<true>(mat->ncol(), rhs, output);
      } else {
        work.template running<false>(mat->ncol(), rhs, output);
      }
    }
  }
  
  template<class Right>
  void adjoint_multiply(const Right& rhs, Workspace& work, Eigen::VectorXd& output) const {
    if (mat->prefer_rows()) {
      if (mat->sparse()) {
        work.template running<true>(mat->nrow(), rhs, output);
      } else {
        work.template running<false>(mat->nrow(), rhs, output);
      }
    } else {
      if (mat->sparse()) {
        work.template direct<true>(mat->ncol(), rhs, output);
      } else {
        work.template direct<false>(mat->ncol(), rhs, output);
      }
    }
  }
  
  Eigen::MatrixXd realize() const {
    Eigen::MatrixXd output(mat->nrow(), mat->ncol());
    tatami::convert_to_dense<false>(mat, output.data());
    return output;
  }
};

//[[Rcpp::export(rng=false)]]
Rcpp::List irlba_tatami(SEXP input, int rank, int nthreads, int seed) {
  auto converted = Rtatami::BoundNumericPointer(input);
  const auto& shared = converted->ptr;
  size_t NR = shared->nrow();
  size_t NC = shared->ncol();
  
  auto trans = tatami::make_DelayedTranspose(std::move(shared));
  TatamiWrapper mat(trans.get(), nthreads);
  
  irlba::EigenThreadScope t(nthreads);
  irlba::Irlba runner;
  runner.set_number(rank);
  runner.set_seed(seed);
  
  Eigen::MatrixXd u;
  Eigen::MatrixXd v;
  Eigen::VectorXd d;
  runner.run(mat, u, v, d);
  
  Rcpp::NumericMatrix out_u(u.rows(), u.cols());
  std::copy(u.data(), u.data() + u.size(), out_u.begin());
  Rcpp::NumericMatrix out_v(v.rows(), v.cols());
  std::copy(v.data(), v.data() + v.size(), out_v.begin());
  Rcpp::NumericVector out_d(d.begin(), d.end());
  
  return Rcpp::List::create(
    Rcpp::Named("u") = out_u,
    Rcpp::Named("v") = out_v,
    Rcpp::Named("d") = out_d
  );
}