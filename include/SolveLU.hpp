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
#include <functional>
#include <iostream>
#include <algorithm>

#include "Exception.hpp"
#include "Matrix.hpp"
#include "LU.hpp"
#include "Pivot.hpp"

#ifndef ANPI_SOLVE_LU_HPP
#define ANPI_SOLVE_LU_HPP

namespace anpi
{
	
namespace fallback 
{

/**
   * Solves a linear equation system with use of the LU decomposition
   *
   * @param[in] A a square matrix 
   * @param[out] x vector of variables to find
   * @param[out] b system's solutions vector.
   *
   * @throws anpi::Exception if a division by zero may occur
   */
template <typename T>
bool solveLU(const anpi::Matrix<T> &A,
             std::vector<T> &x,
             const std::vector<T> &b)
{
    const T eps = std::numeric_limits<T>::epsilon();
    anpi::Matrix<T> LU, L, U;
    std::vector<size_t> permut;
    anpi::lu(A, LU, permut);
    int n = LU.rows();
    anpi::Matrix<T> P = anpi::Matrix<T>(n, n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if ((int)permut[i] == j)
            {
                P[i][j] = 1;
                break;
            }
        }
    }
    anpi::lumpl::unpackDoolittle(LU, L, U);
    std::vector<T> Pb = P * b;
    std::vector<T> y;
    if (abs(L[0][0]) < eps)
    {
        throw anpi::Exception("Can't divide by 0");
    }
    y.push_back(Pb[0] / L[0][0]);
    //Sustitución hacia adelante Ly = Pb y obtener y
    T sum;
    for (int i = 1; i < n; ++i)
    {
        sum = 0;
        if (abs(L[i][i]) < eps)
        {
            throw anpi::Exception("Can't divide by 0");
        }
        for (int j = 0; j < i; ++j)
        {   
            sum += L[i][j] * y[j];
        }
        y.push_back((Pb[i] - sum) / L[i][i]);
    }

    //Sustitución hacia atras Ux = y
    if (abs(U[n - 1][n - 1]) < eps)
    {
        throw anpi::Exception("Can't divide by 0");
    }

    x = std::vector<T>(n, 0);
    x[n - 1] = y[n - 1] / U[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        sum = 0;
        if (abs(U[i][i]) < eps)
        {
            throw anpi::Exception("Can't divide by 0");
        }
        for (int j = i + 1; j < n; ++j)
        {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return true;
}

} // fallback

namespace simd 
{

/**
   * Solves a linear equation system with use of the LU decomposition
   *
   * @param[in] A a square matrix 
   * @param[out] x vector of variables to find
   * @param[out] b system's solutions vector.
   *
   * @throws anpi::Exception if a division by zero may occur
   */
template <typename T>
bool solveLU(const anpi::Matrix<T> &A,
             std::vector<T> &x,
             const std::vector<T> &b)
{
    const T eps = std::numeric_limits<T>::epsilon();
    anpi::Matrix<T> LU, L, U;
    std::vector<size_t> permut;
    anpi::lu(A, LU, permut);
    int n = LU.rows();
    anpi::Matrix<T> P = anpi::Matrix<T>(n, n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if ((int)permut[i] == j)
            {
                P[i][j] = 1;
                break;
            }
        }
    }
    anpi::lumpl::unpackDoolittle(LU, L, U);
    std::vector<T> Pb = P * b;
    std::vector<T> y;
    if (abs(L[0][0]) < eps)
    {
        throw anpi::Exception("Can't divide by 0");
    }
    y.push_back(Pb[0] / L[0][0]);
    //Sustitución hacia adelante Ly = Pb y obtener y
    T sum;
    for (int i = 1; i < n; ++i)
    {
        sum = 0;
        if (abs(L[i][i]) < eps)
        {
            throw anpi::Exception("Can't divide by 0");
        }
        for (int j = 0; j < i; ++j)
        {   
            sum += L[i][j] * y[j];
        }
        y.push_back((Pb[i] - sum) / L[i][i]);
    }

    //Sustitución hacia atras Ux = y
    if (abs(U[n - 1][n - 1]) < eps)
    {
        throw anpi::Exception("Can't divide by 0");
    }

    x = std::vector<T>(n, 0);
    x[n - 1] = y[n - 1] / U[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        sum = 0;
        if (abs(U[i][i]) < eps)
        {
            throw anpi::Exception("Can't divide by 0");
        }
        for (int j = i + 1; j < n; ++j)
        {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return true;
}

} // simd


#ifdef ANPI_ENABLE_SIMD
	namespace simpl = simd;
#else
	namespace simpl = fallback;
#endif

} //anpi

#endif
