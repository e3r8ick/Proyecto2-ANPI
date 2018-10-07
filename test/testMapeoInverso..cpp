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

#include <functional>

#include <cmath>

namespace anpi {
  namespace test {
    
    /// Test the ResistorGrid map
    template<typename T>
    void mapeoInversoTest() {

      // The result
      anpi::indexPair result;

      // teste de indices iguales
      {
        Matrix<float> A = {{1.0,7.0,6.0,4.0},{2.0,17.0,27.0,17.0}};
        std::vector<float> p;
        const size_t idx = 12;


        ResistorGrid a = ResistorGrid(A,p,A);
        result = a.ResistorGrid::indexToNodes(idx);
        BOOST_CHECK(result.row1==1);
        BOOST_CHECK(result.col1==1);
        BOOST_CHECK(result.row2==2);
        BOOST_CHECK(result.col2==1);
      }
    }

  } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( MapeoInverso )

BOOST_AUTO_TEST_CASE(ResistorGrid) 
{
  anpi::test::mapeoInversoTest<float>();
}


BOOST_AUTO_TEST_SUITE_END()
