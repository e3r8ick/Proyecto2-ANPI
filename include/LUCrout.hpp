/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <functional>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_CROUT_HPP
#define ANPI_LU_CROUT_HPP

namespace anpi
{

/**
   * Auxiliary method used to debug LU decomposition.
   *
   * It separates a packed LU matrix into the lower triangular matrix
   * L and the upper triangular matrix U, such that the diagonal of U
   * is composed by 1's.
   */
template <typename T>
void unpackCrout(const Matrix<T> &LU,
                 Matrix<T> &L,
                 Matrix<T> &U)
{

  U = LU;
  int rows = LU.rows();
  int cols = LU.cols();
  if (rows != cols)
  {
    throw anpi::Exception("Matriz debe ser cuadrada.");
  }
  L = Matrix<T>(rows, cols, 0.0);
  U[0][0] = 1;

  for (int i = 0; i < rows; ++i)
  {
    U[i][i] = 1;
    L[i][i] = LU[i][i];
    for (int j = 0; j < i; j++)
    {
      U[i][j] = 0;
      L[i][j] = LU[i][j];
    }
  }
}

/**
   * Decompose the matrix A into a lower triangular matrix L and an
   * upper triangular matrix U.  The matrices L and U are packed into
   * a single matrix LU.  
   *
   * Crout's way of packing assumes a diagonal of
   * 1's in the U matrix.
   *
   * @param[in] A a square matrix 
   * @param[out] LU matrix encoding the L and U matrices
   * @param[out] permut permutation vector, holding the indices of the
   *             original matrix falling into the corresponding element.
   *             For example if permut[5]==3 holds, then the fifth row
   *             of the LU decomposition in fact is dealing with the third
   *             row of the original matrix.
   *
   * @throws anpi::Exception if matrix cannot be decomposed, or input
   *         matrix is not square.
   * @throws anpi::Exception if after permutation there's a 0 in the 
   *         diagonal
   */

template <typename T>
void luCrout(const Matrix<T> &A,
             Matrix<T> &LU,
             std::vector<size_t> &permut)
{
  //pivoteo
  pivot(A, LU, permut);

  //descomposición LU
  int n = A.rows();
  int m = A.cols();
  if (m != n)
  {
    throw anpi::Exception("La matriz debe ser cuadrada");
  }
  if (LU[0][0] == 0)
  {
    throw anpi::Exception("Imposible realizar la descomposición con 0 en la diagonal");
  }
  for (int j = 1; j < n; ++j)
  {
    LU[0][j] = LU[0][j] / LU[0][0];
  }
  T sum;
  for (int j = 1; j < n; ++j)
  {

    for (int i = j; i < n; ++i)
    {
      sum = 0;
      for (int k = 0; k < j; ++k)
      {
        sum += LU[i][k] * LU[k][j];
      }
      LU[i][j] -= sum;
    }

    for (int k = j + 1; k < n; ++k)
    {
      sum = 0;
      for (int i = 0; i < j; ++i)
      {
        sum += LU[j][i] * LU[i][k];
      }
      try
      {
        LU[j][k] = (LU[j][k] - sum) / LU[j][j];
      }
      catch (anpi::Exception& e)
      {
        throw anpi::Exception("Imposible realizar la descomposición con 0 en la diagonal");
      }
    }
      sum = 0;
  }
}

} //anpi

#endif
