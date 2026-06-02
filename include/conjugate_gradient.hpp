#pragma once

#include <sparse_matrix.hpp>

//========================================

template<typename T>
Matrix<T> ConjugateGradientSolve(
	const SparseMatrix<T>& A, 
	const Matrix<T>& b,
	size_t max_iterations = 10
);

//========================================

template<typename T>
Matrix<T> ConjugateGradientSolve(
	const SparseMatrix<T>& A, 
	const Matrix<T>& b,
	size_t max_iterations /*= 10*/
)
{
	Matrix<T> x_prev(b.getHeight(), b.getWidth());
	Matrix<T> r_prev = b;
	Matrix<T> v = r_prev;

	Matrix<T> Av(b.getHeight(), b.getWidth());
	for (size_t i = 0; i < max_iterations; i++)
	{
		auto Av = A * v;
		auto r_prev_sqr = r_prev.dot(r_prev);

		auto t = r_prev_sqr / v.dot(Av);
		auto x = x_prev + v * t;
		auto r = r_prev - Av * t;
		auto s = r.dot(r) / r_prev_sqr;
		v = r + v * s;

		x_prev = x;
		r_prev = r;
	}

	return x_prev;
}

//========================================