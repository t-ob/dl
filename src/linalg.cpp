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
// template <class T>
// Matrix<T> identity_matrix(int n)
// {
//     Matrix m(n, n);
//     for (int i = 0; i < n; ++i)
//         m.set(i, i, 1.0);

//     return m;
// }

template <class T>
Matrix<T> diagonal_matrix(int m, int n, const std::vector<T>& ds)
{
    if (ds.size() > std::min(m, n))
        throw std::domain_error("Too many entries");

    Matrix<T> diag(m, n);
    for (int i = 0; i < ds.size(); ++i)
        diag.set(i, i, ds[i]);

    return diag;
}

template <class T>
Matrix<T> diagonal_matrix(const std::vector<T>& ds) {
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
            double c_ij = 0.0;
            for (int k = 0; k < a.cols(); ++k)
            {
                c_ij += a.get(i, k) * b.get(k, j);
            }
            c.set(i, j, c_ij);
        }
    }

    return c;
}

unsigned long long fib(int n)
{
    Matrix<unsigned long long> f(2, 2);
    f.set(0, 0, 0);
    f.set(0, 1, 1);
    f.set(1, 0, 1);
    f.set(1, 1, 1);

    Matrix<unsigned long long> b(2, 1);
    b.set(0, 0, 0);
    b.set(1, 0, 1);

    while (n > 0) {
        b = f * b;
        n--;
    }

    return b.get(0, 0);
}

int main()
{
    std::vector<int> ds;
    ds.push_back(1);
    ds.push_back(2);
    ds.push_back(3);

    Matrix<int> m = diagonal_matrix(ds);

    repr(m);

    std::cout << fib(99) << std::endl;

    return 0;
}