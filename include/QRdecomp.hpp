#include <cmath>
#include <limits>
#include <iostream>
#include <functional>
#include <vector>
#include <algorithm>
#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_QR_DECOMP_HPP
#define ANPI_QR_DECOMP_HPP

namespace anpi
{

/**
 * calculate the square of a number x
 * */	
template <typename T>
T sqr(T x)
{
	return x * x;
}

/**
 * get the sign of a number x
 * 
*/
template <typename T>
int sgn(T x)
{
	return (T(0) < x) - (x < T(0));
}


/**
 * Perform QR decomposition
 * @param A original matrix
 * @param Q orthogonal matrix
 * @param R right matrix
 * 
*/
template <typename T>
void qr(const anpi::Matrix<T> &A,
		anpi::Matrix<T> &Q,
		anpi::Matrix<T> &R)
{
	int n = A.rows();
	Q = anpi::Matrix<T>(n, n, 0.0);
	R = A;
	std::vector<T> c(n);
	std::vector<T> d(n);
	T scale, sigma, sum, tau;

	for (int k = 0; k < n - 1; ++k)
	{
		scale = 0.0;
		for (int i = k; i < n; ++i)
		{
			scale = std::max(scale, std::abs(R(i, k)));
		}
		if (scale == 0.0)// singular case
		{
			c[k] = d[k] = 0.0;
			throw anpi::Exception("Intentando resolver sistema singular");
		}
		else    //form Qk and Qk *A 
		{
			for (int i = k; i < n; ++i)
			{
				R(i, k) = R(i, k) / scale;
			}
			sum = 0.0;
			for (int i = k; i < n; ++i)
			{
				sum = sum + sqr(R(i, k));
			}
			if (R[k][k] >= 0)
			{
				sigma = sqrt(sum);
			}
			else
			{
				sigma = -sqrt(sum);
			}
			R(k, k) = R(k, k) + sigma;
			c[k] = sigma * R(k, k);
			d[k] = -scale * sigma;
			for (int j = k + 1; j < n; ++j)
			{
				sum = 0.0;
				for (int i = k; i < n; ++i)
				{
					sum = sum + R(i, k) * R(i, j);
				}
				tau = sum / c[k];
				for (int i = k; i < n; i++)
				{
					R(i, j) = R(i, j) - tau * R(i, k);
				}
			}
		}
	}

	d[n - 1] = R(n - 1, n - 1);
	if (d[n - 1] == 0.0)
	{
		throw anpi::Exception("Intentando resolver sistema singular");
	}

	//Form Qt
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			Q(i, j) = 0.0;
		}
		Q(i, i) = 1.0;
	}

	for (int k = 0; k < n - 1; ++k)
	{
		if (c[k] != 0.0)
		{
			for (int j = 0; j < n; ++j)
			{
				sum = 0.0;
				for (int i = k; i < n; i++)
				{
					sum = sum + R(i, k) * Q(i, j);
				}
				sum = sum / c[k];
				for (int i = k; i < n; ++i)
				{
					Q(i, j) = Q(i, j) - sum * R(i, k);
				}
			}
		}
	}
	for (int i = 0; i < n; ++i)
	{
		R(i, i) = d[i];
		for (int j = 0; j < i; ++j)
		{
			R(i, j) = 0.0;
		}
	}

}
/**
 * solve a set of linear aquations Ax=b
 * @param A original matrix
 * @param x vector that we are looking for
 * @param b solution vector
*/
template <typename T>
bool solveQR(const anpi::Matrix<T> &A,
			 std::vector<T> &x,
			 const std::vector<T> &b)
{
	anpi::Matrix<T> Q = anpi::Matrix<T>(A.rows(), A.cols(), 0.0);
	anpi::Matrix<T> R = anpi::Matrix<T>(A.rows(), A.cols(), 0.0);
	anpi::qr(A, Q, R);
	qtMult(b, x, Q);
	rSolve(x, x, R);
	return true;
}


/**
 * perform Qt*b and put the result in x
 * 
*/
template <typename T>
void qtMult(const std::vector<T> &b, std::vector<T> &x, anpi::Matrix<T> &Q)
{
	T sum;	
	int n = Q.rows();
	for (int i = 0; i < n; ++i)
	{
		sum = 0.0;
		for (int j = 0; j < n; ++j)
		{

			sum = sum + Q(i, j) * b[j];
			
		}

		x[i] = sum;
	}

}

/**
 * solve the set of linear equations Rx=b
 * 
*/
template <typename T>
void rSolve(const std::vector<T> &b, std::vector<T> &x, anpi::Matrix<T> &R)
{

	T sum;
	int n = R.rows();
	for (int i = n - 1; i >= 0; --i)
	{
		sum = b[i];
		for (int j = i + 1; j < n; ++j)
		{
			sum = sum - R(i, j) * x[j];
		}
		x[i] = sum / R(i, i);
	}
}
}
#endif