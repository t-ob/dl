#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "linalg.h"

// O(m*n) space
template <class T>
Matrix<T>::Matrix(int r, int c) : _rows(r), _cols(c)
{
    _data.resize(r * c, T());
}

template <class T>
Matrix<T>::Matrix(int r, int c, T x) : _rows(r), _cols(c)
{
    _data.resize(r * c, T());
    for (int i = 0; i < r; ++i)
    {
        _data[i * c + i] = x;
    }
}

// O(1)
template <class T>
void Matrix<T>::set(int i, int j, T x)
{
    _data[_cols * i + j] = x;
}

// O(1)
template <class T>
T Matrix<T>::get(int i, int j) const
{
    return _data[_cols * i + j];
}

// O(1)
template <class T>
int Matrix<T>::rows() const
{
    return _rows;
}

// O(1)
template <class T>
int Matrix<T>::cols() const
{
    return _cols;
}

template <class T>
T Matrix<T>::trace() const
{
    T trace;
    for (int i = 0; i < std::min(_rows, _cols); ++i)
    {
        trace += get(i, i);
    }

    return trace;
}

// O(m*n) space and complexity
template <class T>
Matrix<T> Matrix<T>::transpose() const
{
    Matrix mt(_cols, _rows);
    for (int i = 0; i < _rows; ++i)
    {
        for (int j = 0; j < _cols; ++j)
        {
            mt.set(j, i, get(i, j));
        }
    }

    return mt;
}

// O(n*n) space and O(n) complexity
template <class T>
Matrix<T> identity_matrix(int n)
{
    Matrix<T> m(n, n);
    for (int i = 0; i < n; ++i)
        m.set(i, i, (T)1);

    return m;
}

template <class T>
Matrix<T> diagonal_matrix(int m, int n, const std::vector<T> &ds)
{
    if (ds.size() > std::min(m, n))
        throw std::domain_error("Too many entries");

    Matrix<T> diag(m, n);
    for (int i = 0; i < ds.size(); ++i)
        diag.set(i, i, ds[i]);

    return diag;
}

template <class T>
Matrix<T> diagonal_matrix(const std::vector<T> &ds)
{
    return diagonal_matrix(ds.size(), ds.size(), ds);
}

template <class T>
void repr(Matrix<T> m)
{
    for (int i = 0; i < m.rows(); ++i)
    {
        for (int j = 0; j < m.cols(); ++j)
        {
            std::cout << m.get(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

// O(m*n) space and O(m*n) complexity
template <class T>
Matrix<T> operator+(const Matrix<T> &a, const Matrix<T> &b)
{
    // Matrices be of same shape
    if (a.rows() != b.rows() || a.cols() != b.cols())
        throw std::domain_error("Matrices wrong shape");

    Matrix<T> c(a.rows(), a.cols());
    for (int i = 0; i < a.rows(); ++i)
    {
        for (int j = 0; j < a.cols(); ++j)
        {
            c.set(i, j, a.get(i, j) + b.get(i, j));
        }
    }

    return c;
}

// O(a_n*b_m space) and O(a_n*b_m*a_m) complexity
template <class T>
Matrix<T> operator*(const Matrix<T> &a, const Matrix<T> &b)
{
    // Number of columns in LHS must be equal to rows in RHS
    if (a.cols() != b.rows())
        throw std::domain_error("Incompatible shape");

    Matrix<T> c(a.rows(), b.cols());
    for (int i = 0; i < a.rows(); ++i)
    {
        for (int j = 0; j < b.cols(); ++j)
        {
            T c_ij = (T)0;
            for (int k = 0; k < a.cols(); ++k)
            {
                c_ij += a.get(i, k) * b.get(k, j);
            }
            c.set(i, j, c_ij);
        }
    }

    return c;
}

template <class T>
Matrix<T> pow(const Matrix<T> &m, int exp)
{
    if (m.rows() != m.cols())
        throw std::domain_error("Non-square matrix input");
    if (exp < 1)
        throw std::domain_error("Non-positive exponent");

    Matrix<T> a(m.rows(), m.rows(), (T)1);
    Matrix<T> b = m;
    while (exp > 1)
    {
        if (exp % 2 == 0)
        {
            b = b * b;
            exp /= 2;
        }
        else
        {
            a = a * b;
            exp -= 1;
        }
    }

    return a * b;
}

// Only good up to n = 93 for 64bit ULL!
unsigned long long fib(int n)
{
    Matrix<unsigned long long> f(2, 2);
    f.set(0, 0, 0ULL);
    f.set(0, 1, 1ULL);
    f.set(1, 0, 1ULL);
    f.set(1, 1, 1ULL);

    return pow(f, n).get(0, 1);
}

int main()
{
    for (int i = 1; i < 94; ++i)
        std::cout << i << "\t" << fib(i) << std::endl;

    return 0;
}