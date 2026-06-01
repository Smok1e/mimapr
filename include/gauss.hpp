#include <stdexcept>

#include <matrix.hpp>

//========================================

template<typename T>
Matrix<T> GaussSolve(const Matrix<T>& A, const Matrix<T>& b);

//========================================

template<typename T>
Matrix<T> GaussSolve(const Matrix<T>& A, const Matrix<T>& b)
{
	if (b.getWidth() != 1)
		throw std::invalid_argument("b is not a column-vector");

	if (A.getHeight() != b.getHeight())
		throw std::invalid_argument("coefficient matrix height is different from the column-vector");

	// Прямой ход
	auto n = A.getHeight();

	Matrix<T> Ab = A;
	Ab.resize(n, Ab.getWidth() + b.getWidth());

	for (size_t y = 0; y < n; y++)
		Ab[y][A.getWidth()] = b[y][0];

	for (size_t k = 0; k < n - 1; k++)
	{
		for (size_t i = k + 1; i < n; i++)
		{
			auto coeff = Ab[i][k] / Ab[k][k];

			for (size_t j = k; j < Ab.getWidth(); j++)
				Ab[i][j] -= Ab[k][j] * coeff;
		}
	}

	// Обратный ход
	Matrix<T> solution(Ab.getHeight(), 1);
	for (size_t y = 0; y < n; y++)
		solution[y][0] = Ab[y][A.getWidth()];

	for (int k = n - 1; k >= 0; k--)
	{
		for (int i = n - 1; i > k; i--)
			solution[k][0] -= Ab[k][i] * solution[i][0];

		solution[k][0] /= Ab[k][k];
	}

	return solution;
}

//========================================