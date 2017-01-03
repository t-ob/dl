#ifndef GUARD_linalg_h
#define GUARD_linalg_h

#include <vector>

template <class T>
class Matrix {
    private:
    int _rows;
    int _cols;
    std::vector<T> _data;

    public:
    void set(int, int, T);
    T get(int, int) const;
    int rows() const;
    int cols() const;

    T trace() const;
    Matrix transpose() const;

    Matrix(int, int);
    Matrix(int, int, T);
};  

#endif
