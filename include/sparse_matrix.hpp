#pragma once

#include <stdexcept>
#include <vector>

#include <matrix.hpp>

//========================================

template<typename T>
class SparseMatrix
{
public:
	SparseMatrix(const Matrix<T>& mat, double eps = 1e-9);
	SparseMatrix(const SparseMatrix<T>& copy) = delete;

	Matrix<T> operator*(const Matrix<T>& vec) const;

private:
	size_t m_width;
	size_t m_height;
	std::vector<T> m_values;
	std::vector<T> m_col_index;
	std::vector<T> m_row_index;

};

//========================================

template<typename T>
SparseMatrix<T>::SparseMatrix(const Matrix<T>& mat, double eps /*= 1e-9*/):
	m_width(mat.getWidth()),
	m_height(mat.getHeight())
{
	m_row_index.push_back(0);

	for (size_t row = 0; row < m_height; row++)
	{
		size_t nonzero_count = 0;

		for (size_t col = 0; col < m_width; col++)
		{
			if (abs(mat[row][col]) > eps)
			{
				nonzero_count++;
				m_values.push_back(mat[row][col]);
				m_col_index.push_back(col);
			}
		}

		m_row_index.push_back(m_row_index.back() + nonzero_count);
	}
}

//========================================

template<typename T>
Matrix<T> SparseMatrix<T>::operator*(const Matrix<T>& vec) const
{
	if (m_width != vec.getHeight() || m_height != vec.getHeight())
		throw std::invalid_argument("incompatiable vector");

	if (vec.getWidth() != 1)
		throw std::invalid_argument("right hand side operand is not a vector");

	Matrix<T> product(vec.getHeight(), vec.getWidth());
	for (size_t row = 0; row < m_height; row++)
	{
		size_t cols_begin = m_row_index[row    ];
		size_t cols_end   = m_row_index[row + 1];
		size_t col_count = cols_end - cols_begin;

		for (size_t i = 0; i < col_count; i++)
			product[row][0] += m_values[cols_begin + i] * vec[m_col_index[cols_begin + i]][0];
	}

	return product;
}

//========================================