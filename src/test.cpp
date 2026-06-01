#include <iostream>
#include <print>

#include <matrix.hpp>
#include <cholesky.hpp>

//========================================

int main()
{
	Matrix<double> A({
		{ 25, 15, -5 },
		{ 15, 18,  0 },
		{ -5,  0, 11 }
	});

	Matrix<double> b({
		{ 1 },
		{ 123 },
		{ 2 }
	});

	auto L = CholeskyDecompose(A);
	auto x = CholeskySolve(L, b);

	std::cout << x << std::endl;

	for (size_t i = 0; i < A.getHeight(); i++)
	{
		double sum = 0;
		for (size_t j = 0; j < A.getWidth(); j++)
			sum += A[i][j] * x[j][0];

		std::println("{}", sum);
	}
}

//========================================