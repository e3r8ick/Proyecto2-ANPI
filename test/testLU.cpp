/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <boost/test/unit_test.hpp>

#include "LUCrout.hpp"
#include "LUDoolittle.hpp"
#include "Invert.hpp"
#include "SolveLU.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>

namespace anpi {
  namespace test {
    
    /// Test the given closed root finder
    template<typename T>
    void luTest(const std::function<void(const Matrix<T>&,
                                         Matrix<T>&,
                                         std::vector<size_t>&)>& decomp,
                const std::function<void(const Matrix<T>&,
                                         Matrix<T>&,
                                         Matrix<T>&)>& unpack) {

      // The result
      Matrix<T> LU;

      // Test if a non-square matrix is successfully detected
      {
        Matrix<T> A = {{1,7,6,4},{2,17,27,17}};
        std::vector<size_t> p;
        try {
          decomp(A,LU,p);
          BOOST_CHECK_MESSAGE(false,"Rectangular matrix not properly catched");
        }
        catch(anpi::Exception& exc) {
          BOOST_CHECK_MESSAGE(true,"Rectangular matrix properly detected");
        }
      }

      // Test pivoting
      {
        anpi::Matrix<T> A = { {-1,-2,1,2},{ 2, 0,1,2},{-1,-1,0,1},{ 1, 1,1,1} };
        std::vector<size_t> p;
        decomp(A,LU,p);

        std::vector<size_t> gp= {1,0,3,2};
        BOOST_CHECK(gp==p);
      }
      
      // Test decomposition
      {
        // same matrix as before, but already permuted to force a clean decomposition
        anpi::Matrix<T> A = { { 2, 0,1,2},{-1,-2,1,2},{ 1, 1,1,1},{-1,-1,0,1} };
        std::vector<size_t> p;
        decomp(A,LU,p);
        Matrix<T> L,U;
        unpack(LU,L,U);
        Matrix<T> Ar=L*U;

        const T eps = std::numeric_limits<T>::epsilon();
        BOOST_CHECK(Ar.rows()==A.rows());
        BOOST_CHECK(Ar.cols()==A.cols());
        
        for (size_t i=0;i<Ar.rows();++i) {
          for (size_t j=0;j<Ar.cols();++j) {
            BOOST_CHECK(std::abs(Ar(i,j)-A(i,j)) < eps);
          }
        }
      }

      //test invert
      {
        anpi::Matrix<T> A = { { 2, 0,1,2},{-1,-2,1,2},{ 1, 1,1,1},{-1,-1,0,1} };
        anpi::Matrix<T> Ai;
        invert(A, Ai);

        const T eps = std::numeric_limits<T>::epsilon();
        anpi::Matrix<T> I = A*Ai;

        
        for (int i = 0; i < (int) A.rows(); ++i){
          for (int j = 0; j < (int) A.cols(); ++j){
            if (i == j){
              BOOST_CHECK(abs(I[i][j] - 1 ) < eps);
            }else{
              BOOST_CHECK(abs(I[i][j]) < eps);
            }
          }
        }
        
      }

      //Test solve LU
      {
        anpi::Matrix<float> A = {{-1, -2, 1},
                           {2, 0, 1},
                           {-1, -1, 0}};
        std::vector<float> b = {1,0,-1};
        std::vector<float> x;
        simpl::solveLU(A, x, b);

        const T eps = std::numeric_limits<T>::epsilon();
        //primera solución -3
        BOOST_CHECK(abs(x[0]+3) < eps);
        //primera solución 4
        BOOST_CHECK(abs(x[1]-4) < eps);
        //primera solución 6
        BOOST_CHECK(abs(x[2]-6) < eps);
      }

      //Test solve LU with an error 
      {
        anpi::Matrix<float> A = {{0, 0, 0},
                           {0, 0, 0},
                           {0, 0, 0}};
        std::vector<float> b = {1,0,-1};
        std::vector<float> x;

        try {
          simpl::solveLU(A, x, b);
          BOOST_CHECK_MESSAGE(false,"Divide by 0 not catch");
        }
        catch(anpi::Exception& exc) {
          BOOST_CHECK_MESSAGE(true,"Divide by 0 catch");
        }
      }
    }

  } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( LU )

BOOST_AUTO_TEST_CASE(Doolittle) 
{
  anpi::test::luTest<float>(anpi::lumpl::luDoolittle<float>,
                            anpi::lumpl::unpackDoolittle<float>);
  anpi::test::luTest<double>(anpi::lumpl::luDoolittle<double>,
                             anpi::lumpl::unpackDoolittle<double>);
}
/*
BOOST_AUTO_TEST_CASE(Crout) 
{
  anpi::test::luTest<float>(anpi::luCrout<float>,anpi::unpackCrout<float>);
  anpi::test::luTest<double>(anpi::luCrout<double>,anpi::unpackCrout<double>);
}
*/

BOOST_AUTO_TEST_SUITE_END()
