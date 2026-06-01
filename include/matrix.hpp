#pragma once

#include <ostream>
#include <print>
#include <ranges>

//========================================

template<typename T = double>
class Matrix
{
public:
	Matrix(size_t height, size_t width);

	template<size_t Height, size_t Width>
	explicit Matrix(const T (&data)[Height][Width]);

	Matrix(const Matrix<T>& copy);
	Matrix(Matrix<T>&& move) noexcept = default;
	~Matrix();

	size_t getWidth() const;
	size_t getHeight() const;

	// Get row pointer
	T* get(size_t y);

	// Get const row pointer
	const T* get(size_t y) const;

	// Get element reference
	T& get(size_t y, size_t x);

	// Get const element reference
	const T& get(size_t y, size_t x) const;

	// Get row pointer
	T* operator()(size_t y);

	// Get const row pointer
	const T* operator()(size_t y) const;

	// Get row pointer
	T* operator[](size_t y);

	// Get const row pointer
	const T* operator[](size_t y) const;

	// Get element reference
	T& operator()(size_t y, size_t x);

	// Get const element reference
	const T& operator()(size_t y, size_t x) const;

	T* begin();
	const T* begin() const;

	T* end();
	const T* end() const;

	void resize(size_t new_height, size_t new_width);

private:
	size_t m_width;
	size_t m_height;

	T* m_data = nullptr;

};

//========================================

template<typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix);

//========================================

template<typename T>
Matrix<T>::Matrix(size_t height, size_t width):
	m_width(width),
	m_height(height)
{
	m_data = new T[width * height](T {});
}

template<typename T>
template<size_t Height, size_t Width>
Matrix<T>::Matrix(const T (&data)[Height][Width]):
	Matrix(Height, Width)
{
	for (size_t y = 0; y < Height; y++)
		std::copy(data[y], data[y] + Width, m_data + y * Width);
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& copy):
	Matrix(copy.m_height, copy.m_width)
{
	std::copy(copy.begin(), copy.end(), begin());
}

template<typename T>
Matrix<T>::~Matrix()
{
	delete[] m_data;
}

//========================================

template<typename T>
size_t Matrix<T>::getWidth() const
{
	return m_width;
}

template<typename T>
size_t Matrix<T>::getHeight() const
{
	return m_height;
}

template<typename T>
T* Matrix<T>::get(size_t y)
{
	return m_data + y * m_width;
}

template<typename T>
const T* Matrix<T>::get(size_t y) const
{
	return m_data + y * m_width;
}

template<typename T>
T& Matrix<T>::get(size_t y, size_t x)
{
	return get(y)[x];
}

template<typename T>
const T& Matrix<T>::get(size_t y, size_t x) const
{
	return get(y)[x];
}

template<typename T>
T& Matrix<T>::operator()(size_t y, size_t x)
{
	return get(y, x);
}

template<typename T>
const T& Matrix<T>::operator()(size_t y, size_t x) const
{
	return get(y, x);
}

template<typename T>
T* Matrix<T>::operator()(size_t y)
{
	return get(y);
}

template<typename T>
const T* Matrix<T>::operator()(size_t y) const
{
	return get(y);
}

template<typename T>
T* Matrix<T>::operator[](size_t y)
{
	return get(y);
}

template<typename T>
const T* Matrix<T>::operator[](size_t y) const
{
	return get(y);
}

template<typename T>
T* Matrix<T>::begin()
{
	return m_data;
}

template<typename T>
const T* Matrix<T>::begin() const
{
	return m_data;
}

template<typename T>
T* Matrix<T>::end()
{
	return m_data + m_width * m_height;
}

template<typename T>
const T* Matrix<T>::end() const
{
	return m_data + m_width * m_height;
}

//========================================

template<typename T>
void Matrix<T>::resize(size_t new_height, size_t new_width)
{
	T* new_data = new T[new_width * new_height](static_cast<T>(0));

	auto min_width  = std::min(m_width,  new_width );
	auto min_height = std::min(m_height, new_height);

	for (size_t y = 0; y < min_height; y++)
		std::copy(get(y), get(y) + min_width, new_data + y * new_width);

	std::swap(m_data, new_data);
	std::swap(m_width, new_width);
	std::swap(m_height, new_height);
	delete[] new_data;
}

//========================================

template<typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix)
{
	for (size_t y = 0; y < matrix.getHeight(); y++)
	{
		for (size_t x = 0; x < matrix.getWidth(); x++)
			std::print(stream, "{: 5.2} ", static_cast<double>(matrix(y, x)));

		std::println(stream, "");
	}

	return stream;
}

//========================================