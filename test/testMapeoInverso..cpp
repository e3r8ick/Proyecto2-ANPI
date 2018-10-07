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
      std::size_t result;

      // teste de indices iguales
      {
        Matrix<float> A = {{1.0,7.0,6.0,4.0},{2.0,17.0,27.0,17.0}};
        std::vector<float> p;
        const size_t idx = 1;


        ResistorGrid a = ResistorGrid(A,p,A);
        a.ResistorGrid::indexToNodes(idx);
        BOOST_CHECK();
      }
    }

  } // test
}  // anpi

BOOST_AUTO_TEST_SUITE( Mapeo )

BOOST_AUTO_TEST_CASE(ResistorGrid) 
{
  anpi::test::mapeoInversoTest<float>();
}


BOOST_AUTO_TEST_SUITE_END()
