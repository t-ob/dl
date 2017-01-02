#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "linalg.h"

// O(m*n) space
Matrix::Matrix(int r, int c) : _rows(r), _cols(c)
{
    _data.resize(r * c, 0.0);
}

// O(1)
void Matrix::set(int i, int j, double x)
{
    _data[_cols * i + j] = x;
}

// O(1)
double Matrix::get(int i, int j) const
{
    return _data[_cols * i + j];
}

// O(1)
int Matrix::rows() const
{
    return _rows;
}

// O(1)
int Matrix::cols() const
{
    return _cols;
}

// O(m*n) space and complexity
Matrix Matrix::transpose() const
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
Matrix identity_matrix(int n)
{
    Matrix m(n, n);
    for (int i = 0; i < n; ++i)
        m.set(i, i, 1.0);

    return m;
}

Matrix diagonal_matrix(int m, int n, const std::vector<double>& ds)
{
    if (ds.size() > std::min(m, n))
        throw std::domain_error("Too many entries");

    Matrix diag(m, n);
    for (int i = 0; i < ds.size(); ++i)
        diag.set(i, i, ds[i]);

    return diag;
}

Matrix diagonal_matrix(const std::vector<double>& ds) {
    return diagonal_matrix(ds.size(), ds.size(), ds);
}

void repr(Matrix m)
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
Matrix operator+(const Matrix &a, const Matrix &b)
{
    // Matrices be of same shape
    if (a.rows() != b.rows() || a.cols() != b.cols())
        throw std::domain_error("Matrices wrong shape");

    Matrix c(a.rows(), a.cols());
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
Matrix operator*(const Matrix &a, const Matrix &b)
{
    // Number of columns in LHS must be equal to rows in RHS
    if (a.cols() != b.rows())
        throw std::domain_error("Incompatible shape");

    Matrix c(a.rows(), b.cols());
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

double fib(int n)
{
    Matrix f(2, 2);
    f.set(0, 0, 0.0);
    f.set(0, 1, 1.0);
    f.set(1, 0, 1.0);
    f.set(1, 1, 1.0);

    Matrix b(2, 1);
    b.set(0, 0, 0.0);
    b.set(1, 0, 1.0);

    while (n > 0) {
        b = f * b;
        n--;
    }

    return b.get(0, 0);
}

int main()
{
    std::vector<double> ds;
    ds.push_back(1.0);
    ds.push_back(2.0);
    ds.push_back(3.0);

    Matrix m = diagonal_matrix(ds);

    repr(m);

    return 0;
}