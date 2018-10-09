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

#include "ResistorGrid.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <string>
#include <functional>

#include <cmath>

namespace anpi {
  namespace test {
    
    /// Test the ResistorGrid map
    template<typename T>
    void mapeoTest() {

      // The result
      std::size_t result,result3;

      // test mapeo y mapeo inverso
      {
        Matrix<float> A = {{1.0,7.0,6.0,4.0},{2.0,17.0,27.0,17.0},{1.0,7.0,6.0,4.0},{2.0,17.0,27.0,17.0}};
        std::vector<float> p;
        const size_t row1 = 1;
        const size_t col1 = 1;
        const size_t row2 = 2;
        const size_t col2 = 1;


        ResistorGrid a = ResistorGrid(A,p,A);
        result = a.ResistorGrid::nodesToIndex(row1, col1, row2, col2);
        anpi::indexPair result2;
        result2 = a.ResistorGrid::indexToNodes(result);
        BOOST_CHECK(result2.row1==row1);
        BOOST_CHECK(result2.col1==col1);
        BOOST_CHECK(result2.row2==row2);
        BOOST_CHECK(result2.col2==col2);
        result3 = a.ResistorGrid::nodesToIndex(result2.row1, result2.col1, result2.row2, result2.col2);
        BOOST_CHECK(result3==result);
      }

      // teste de indices iguales
      {
        Matrix<float> A = {{1.0,7.0,6.0,4.0},{2.0,17.0,27.0,17.0}};
        std::vector<float> p;
        const size_t col1 = 1;
        const size_t col2 = 1;
        const size_t row1 = 1;
        const size_t row2 = 1;


        ResistorGrid a = ResistorGrid(A,p,A);
        try{
          a.ResistorGrid::nodesToIndex(row1, col1, row2, col2);
          BOOST_CHECK_MESSAGE(false,"Sames indices not properly catched");
        }
         catch(anpi::Exception& exc) {
          BOOST_CHECK_MESSAGE(true,"Same indices properly catched");
        }
      }

      // teste de indices no adyacentes
      {
        Matrix<float> A = {{1.0,7.0,6.0,4.0},{2.0,17.0,27.0,17.0}};
        std::vector<float> p;
        const size_t col1 = 1;
        const size_t col2 = 1;
        const size_t row1 = 1;
        const size_t row2 = 5;


        ResistorGrid a = ResistorGrid(A,p,A);
        try{
          a.ResistorGrid::nodesToIndex(row1, col1, row2, col2);
          BOOST_CHECK_MESSAGE(false,"Not adjacent indices adyacentes not properly catched");
        }
         catch(anpi::Exception& exc) {
          BOOST_CHECK_MESSAGE(true,"Not adjacent indices properly catched");
        }
      }

      // teste de indices negativos
      {
        Matrix<float> A = {{1.0,7.0,6.0,4.0},{2.0,17.0,27.0,17.0}};
        std::vector<float> p;
        const size_t col1 = -1;
        const size_t col2 = 1;
        const size_t row1 = 1;
        const size_t row2 = 2;


        ResistorGrid a = ResistorGrid(A,p,A);
        try{
          a.ResistorGrid::nodesToIndex(row1, col1, row2, col2);
          BOOST_CHECK_MESSAGE(false,"Negatives indices not properly catched");
        }
         catch(anpi::Exception& exc) {
          BOOST_CHECK_MESSAGE(true,"Negatives indices properly catched");
        }
      }


      //test de carga de imagen
      {
        Matrix<float> A = {{1.0,7.0,6.0,4.0},{2.0,17.0,27.0,17.0}};
        std::vector<float> p;
        ResistorGrid a = ResistorGrid(A,p,A);

        
        a.ResistorGrid::build("/mapa4x5.png");
        BOOST_CHECK(a.getRawMap().rows() == 4);
        BOOST_CHECK(a.getRawMap().cols() == 5);
      }


    }

  } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( Mapeo )

BOOST_AUTO_TEST_CASE(ResistorGrid) 
{
  anpi::test::mapeoTest<float>();
}


BOOST_AUTO_TEST_SUITE_END()
