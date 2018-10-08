/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Pablo Alvarado
 * @date   29.12.2017
 */


#include <boost/test/unit_test.hpp>


#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <ctime>

/**
 * Unit tests for the LU class
 */
#include "benchmarkFramework.hpp"
#include "Matrix.hpp"
#include "Allocator.hpp"
#include "LUDoolittle.hpp"
#include "LUCrout.hpp"

BOOST_AUTO_TEST_SUITE( MatrixLU )

/// Benchmark for addition operations
template<typename T>
class benchLU {
protected:
  /// Maximum allowed size for the square matrices
  const size_t _maxSize;

  /// A large matrix holding 
  anpi::Matrix<T> _data;

  /// State of the benchmarked evaluation
  anpi::Matrix<T> LU;
  anpi::Matrix<T> A = {{1,7,6,4},{2,17,27,17}};
  std::vector<size_t> p;
public:
  /// Construct
  benchLU(const size_t maxSize)
    : _maxSize(maxSize),_data(maxSize,maxSize,anpi::DoNotInitialize) {

    size_t idx=0;
    for (size_t r=0;r<_maxSize;++r) {
      for (size_t c=0;c<_maxSize;++c) {
        _data(r,c)=idx++;
      }
    }
  }

  /// Prepare the evaluation of given size
  void prepare(const size_t size) {
    assert (size<=this->_maxSize);
    this->A=std::move(anpi::Matrix<T>(size,size,_data.data()));
    this->LU=this->A;
  }
};

/// Provide the evaluation Doolittle method (float)
template<typename T>
class benchDoolittleFloat : public benchLU<T> {
public:
  /// Constructor
  benchDoolittleFloat(const size_t n) : benchLU<T>(n) { }
  
  // Evaluate Doolittle
  inline void eval() {
    anpi::fallback::luDoolittle(this->A,this->LU, this->p);
  }
};

/// Provide the evaluation Doolittle method (float) with SIMD optimization
template<typename T>
class benchDoolittleFloatSIMD : public benchLU<T> {
public:
  /// Constructor
  benchDoolittleFloatSIMD(const size_t n) : benchLU<T>(n) { }
  
  // Evaluate Doolittle
  inline void eval() {
    anpi::simd::luDoolittle(this->A,this->LU, this->p);
  }
};

/// Provide the evaluation Doolittle method (Double)
template<typename T>
class benchDoolittleDouble : public benchLU<T> {
public:
  /// Constructor
  benchDoolittleDouble(const size_t n) : benchLU<T>(n) { }
  
  // Evaluate add in-place
  inline void eval() {
    anpi::fallback::luDoolittle(this->A,this->LU, this->p);
  }
};

/// Provide the evaluation Doolittle method (Double) with SIMD optimization
template<typename T>
class benchDoolittleDoubleSIMD : public benchLU<T> {
public:
  /// Constructor
  benchDoolittleDoubleSIMD(const size_t n) : benchLU<T>(n) { }
  
  // Evaluate add in-place
  inline void eval() {
    anpi::simd::luDoolittle(this->A,this->LU, this->p);
  }
};


/**
 * Instantiate and test the methods of the Matrix class
 */
BOOST_AUTO_TEST_CASE( LU ) {

  std::vector<size_t> sizes = {  24,  32,  48,  64,
                                 96, 128, 192, 256};//para que dure menos tiempos se quitaron elementos

  const size_t n=sizes.back();
  const size_t repetitions=100;
  std::vector<anpi::benchmark::measurement> times;
/*
  {
    benchCroutFloat<float>  bc(n);

    // Measure Crout float
    ANPI_BENCHMARK(sizes,repetitions,times,bc);
    
    ::anpi::benchmark::write("Crout_float.txt",times);
    ::anpi::benchmark::plotRange(times,"Crout (float)","r");
  }

  {
    benchCroutDouble<double>  bc(n);

    // Measure Crout double
    ANPI_BENCHMARK(sizes,repetitions,times,bc);
    
    ::anpi::benchmark::write("Crout_double.txt",times);
    ::anpi::benchmark::plotRange(times,"Crout (double) ","g");
  }
*/
  {
    benchDoolittleFloat<float> bd(n);

    // Measure Doolittle float
    ANPI_BENCHMARK(sizes,repetitions,times,bd);

    ::anpi::benchmark::write("Doolittle_float.txt",times);
    ::anpi::benchmark::plotRange(times,"Doolittle (float)","b");
  }

  {
    benchDoolittleDouble<double> bd(n);

    // Measure Doolittle double
    ANPI_BENCHMARK(sizes,repetitions,times,bd);

    ::anpi::benchmark::write("Doolittle_double.txt",times);
    ::anpi::benchmark::plotRange(times,"Doolittle (double)","m");
  }
  
    {
    benchDoolittleDoubleSIMD<float> bd(n);

    // Measure Doolittle float
    ANPI_BENCHMARK(sizes,repetitions,times,bd);

    ::anpi::benchmark::write("Doolittle_float_simd.txt",times);
    ::anpi::benchmark::plotRange(times,"Doolittle SIMD (float)","b");
  }

  {
    benchDoolittleDoubleSIMD<double> bd(n);

    // Measure Doolittle double
    ANPI_BENCHMARK(sizes,repetitions,times,bd);

    ::anpi::benchmark::write("Doolittle_double_simd.txt",times);
    ::anpi::benchmark::plotRange(times,"Doolittle SIMD (double)","m");
  }
  
#if 0
  
  {
    benchCrout<double>  bf(n);

    // Measure  Crout
    ANPI_BENCHMARK(sizes,repetitions,times,bf);
    
    ::anpi::benchmark::write("crout_on_copy_double.txt",times);
    ::anpi::benchmark::plotRange(times,"On-copy (double)","g");
  }

  {
    benchDoolittle<double> bd(n);

    // Measure Doolittle
    ANPI_BENCHMARK(sizes,repetitions,times,bd);

    ::anpi::benchmark::write("Doolittle_double.txt",times);
    ::anpi::benchmark::plotRange(times,"Doolittle (double)","m");
  }
#endif
  
  ::anpi::benchmark::show();
}
  
BOOST_AUTO_TEST_SUITE_END()
