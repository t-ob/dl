#ifndef GUARD_linalg_h
#define GUARD_linalg_h

#include <vector>

class Matrix {
    private:
    int _rows;
    int _cols;
    std::vector<double> _data;

    public:
    void set(int, int, double);
    double get(int, int) const;
    int rows() const;
    int cols() const;

    Matrix transpose() const;

    Matrix(int, int);
};  

#endif
