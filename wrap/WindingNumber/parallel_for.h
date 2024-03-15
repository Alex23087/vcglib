// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PARALLEL_FOR_H
#define IGL_PARALLEL_FOR_H
#  define IGL_INLINE inline
#include <functional>

//#warning "Defining IGL_PARALLEL_FOR_FORCE_SERIAL"
//#define IGL_PARALLEL_FOR_FORCE_SERIAL

namespace igl
{
  /// Functional implementation of a basic, open-mp style, parallel
  /// for loop. If the inner block of a for-loop can be rewritten/encapsulated in
  /// a single (anonymous/lambda) function call `func` so that the serial code
  /// looks like:
  ///
  /// \code{cpp}
  ///     for(int i = 0;i<loop_size;i++)
  ///     {
  ///       func(i);
  ///     }
  /// \endcode
  ///
  /// then `parallel_for(loop_size,func,min_parallel)` will use as many threads as
  /// available on the current hardware to parallelize this for loop so long as
  /// loop_size<min_parallel, otherwise it will just use a serial for loop.
  ///
  /// Often if your code looks like:
  ///
  /// \code{cpp}
  ///     for(int i = 0;i<loop_size;i++)
  ///     {
  ///       …
  ///     }
  /// \endcode
  ///
  /// Then you can make a minimal two-line change to parallelize it:
  ///
  /// \code{cpp}
  ///     //for(int i = 0;i<loop_size;i++)
  ///     parallel_for(loop_size,[&](int i)
  ///     {
  ///       …
  ///     }
  ///     ,1000);
  /// \endcode
  ///
  /// @param[in] loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
  /// @param[in] func  function handle taking iteration index as only argument to compute
  ///     inner block of for loop I.e. for(int i ...){ func(i); }
  /// @param[in] min_parallel  min size of loop_size such that parallel (non-serial)
  ///     thread pooling should be attempted {0}
  /// @return true iff thread pool was invoked
  template<typename Index, typename FunctionType >
  inline bool parallel_for(
    const Index loop_size,
    const FunctionType & func,
    const size_t min_parallel=0);
  /// Functional implementation of an open-mp style, parallel for loop with
  /// accumulation. For example, serial code separated into n chunks (each to be
  /// parallelized with a thread) might look like:
  ///
  /// \code{cpp}
  ///     Eigen::VectorXd S;
  ///     const auto & prep_func = [&S](int n){ S = Eigen:VectorXd::Zero(n); };
  ///     const auto & func = [&X,&S](int i, int t){ S(t) += X(i); };
  ///     const auto & accum_func = [&S,&sum](int t){ sum += S(t); };
  ///     prep_func(n);
  ///     for(int i = 0;i<loop_size;i++)
  ///     {
  ///       func(i,i%n);
  ///     }
  ///     double sum = 0;
  ///     for(int t = 0;t<n;t++)
  ///     {
  ///       accum_func(t);
  ///     }
  /// \endcode
  ///
  /// @param[in] loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
  /// @param[in] prep_func function handle taking n >= number of threads as only
  ///     argument
  /// @param[in] func  function handle taking iteration index i and thread id t as only
  ///     arguments to compute inner block of for loop I.e.
  ///     for(int i ...){ func(i,t); }
  /// @param[in] accum_func  function handle taking thread index as only argument, to be
  ///     called after all calls of func, e.g., for serial accumulation across
  ///     all n (potential) threads, see n in description of prep_func.
  /// @param[in] min_parallel  min size of loop_size such that parallel (non-serial)
  ///     thread pooling should be attempted {0}
  /// @return true iff thread pool was invoked
  template<
    typename Index,
    typename PrepFunctionType,
    typename FunctionType,
    typename AccumFunctionType
    >
  inline bool parallel_for(
    const Index loop_size,
    const PrepFunctionType & prep_func,
    const FunctionType & func,
    const AccumFunctionType & accum_func,
    const size_t min_parallel=0);


// Implementation
inline unsigned int default_num_threads(unsigned int user_num_threads=0) {
    // Thread-safe initialization using Meyers' singleton
    class MySingleton {
    public:
        static MySingleton &instance(unsigned int force_num_threads) {
            static MySingleton instance(force_num_threads);
            return instance;
        }
        
        unsigned int get_num_threads() const { return m_num_threads; }
        
    private:
        static const char* getenv_nowarning(const char* env_var)
        {
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
            return std::getenv(env_var);
#ifdef _MSC_VER
#pragma warning(pop)
#endif
        }
        
        MySingleton(unsigned int force_num_threads) {
            // User-defined default
            if (force_num_threads) {
                m_num_threads = force_num_threads;
                return;
            }
            // Set from env var
            if (const char *env_str = getenv_nowarning("IGL_NUM_THREADS")) {
                const int env_num_thread = atoi(env_str);
                if (env_num_thread > 0) {
                    m_num_threads = static_cast<unsigned int>(env_num_thread);
                    return;
                }
            }
            // Guess from hardware
            const unsigned int hw_num_threads = std::thread::hardware_concurrency();
            if (hw_num_threads) {
                m_num_threads = hw_num_threads;
                return;
            }
            // Fallback when std::thread::hardware_concurrency doesn't work
            m_num_threads = 8u;
        }
        
        unsigned int m_num_threads = 0;
    };
    
    return MySingleton::instance(user_num_threads).get_num_threads();
}
} // end namespace igl
#include <cmath>
#include <cassert>
#include <thread>
#include <vector>
#include <algorithm>

template<typename Index, typename FunctionType >
inline bool igl::parallel_for(
  const Index loop_size,
  const FunctionType & func,
  const size_t min_parallel)
{
  using namespace std;
  // no op preparation/accumulation
  const auto & no_op = [](const size_t /*n/t*/){};
  // two-parameter wrapper ignoring thread id
  const auto & wrapper = [&func](Index i,size_t /*t*/){ func(i); };
  return parallel_for(loop_size,no_op,wrapper,no_op,min_parallel);
}

template<
  typename Index,
  typename PreFunctionType,
  typename FunctionType,
  typename AccumFunctionType>
inline bool igl::parallel_for(
  const Index loop_size,
  const PreFunctionType & prep_func,
  const FunctionType & func,
  const AccumFunctionType & accum_func,
  const size_t min_parallel)
{
  assert(loop_size>=0);
  if(loop_size==0) return false;
  // Estimate number of threads in the pool
  // http://ideone.com/Z7zldb
#ifdef IGL_PARALLEL_FOR_FORCE_SERIAL
  const size_t nthreads = 1;
#else
  const size_t nthreads = igl::default_num_threads();
#endif
  if(loop_size<min_parallel || nthreads<=1)
  {
    // serial
    prep_func(1);
    for(Index i = 0;i<loop_size;i++) func(i,0);
    accum_func(0);
    return false;
  }else
  {
    // Size of a slice for the range functions
    Index slice =
      std::max(
        (Index)std::round((loop_size+1)/static_cast<double>(nthreads)),(Index)1);

    // [Helper] Inner loop
    const auto & range = [&func](const Index k1, const Index k2, const size_t t)
    {
      for(Index k = k1; k < k2; k++) func(k,t);
    };
    prep_func(nthreads);
    // Create pool and launch jobs
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    // Inner range extents
    Index i1 = 0;
    Index i2 = std::min(0 + slice, loop_size);
    {
      size_t t = 0;
      for (; t+1 < nthreads && i1 < loop_size; ++t)
      {
        pool.emplace_back(range, i1, i2, t);
        i1 = i2;
        i2 = std::min(i2 + slice, loop_size);
      }
      if (i1 < loop_size)
      {
        pool.emplace_back(range, i1, loop_size, t);
      }
    }
    // Wait for jobs to finish
    for (std::thread &t : pool) if (t.joinable()) t.join();
    // Accumulate across threads
    for(size_t t = 0;t<nthreads;t++)
    {
      accum_func(t);
    }
    return true;
  }
}

//#ifndef IGL_STATIC_LIBRARY
//#include "parallel_for.cpp"
//#endif
#endif