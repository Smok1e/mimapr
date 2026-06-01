#include <stdexcept>
#include <cmath>

#include <matrix.hpp>

//========================================

template<typename T>
Matrix<T> CholeskyDecompose(const Matrix<T>& A);

template<typename T>
Matrix<T> CholeskySolve(const Matrix<T>& L, const Matrix<T>& b);

//========================================

template<typename T>
Matrix<T> CholeskyDecompose(const Matrix<T>& A)
{
	if (A.getWidth() != A.getHeight())
		throw std::invalid_argument("can't decompose a non-square matrix");

	auto N = A.getWidth();

	Matrix<T> L(N, N);
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j <= i; j++)
		{
			T sum {};

			if (i == j)
			{
				for (size_t k = 0; k < j; k++)
					sum += L[j][k] * L[j][k];

				L[j][j] = std::sqrt(A[j][j] - sum);
			}

			else
			{
				for (size_t k = 0; k < j; k++)
					sum += L[i][k] * L[j][k];

				L[i][j] = (A[i][j] - sum) / L[j][j];
			}
		}
	}

	return L;
}

template<typename T>
Matrix<T> CholeskySolve(const Matrix<T>& L, const Matrix<T>& b)
{
	if (b.getWidth() != 1)
		throw std::invalid_argument("b is not a column-vector");

	if (L.getHeight() != b.getHeight())
		throw std::invalid_argument("coefficient matrix height is different from the column-vector");

	auto N = b.getHeight();

	// Ly = b
	Matrix<T> y = b;
	for (size_t k = 0; k < N; k++)
	{
		for (size_t i = 0; i < k; i++)
			y[k][0] -= L[k][i] * y[i][0];

		y[k][0] /= L[k][k];
	}

	auto& x = y;

	// L^Tx = y
	for (int k = N - 1; k >= 0; k--)
	{
		for (int i = N - 1; i >= k + 1; i--)
			x[k][0] -= L[i][k] * x[i][0];

		x[k][0] /= L[k][k];
	}

	return x;
}

//========================================