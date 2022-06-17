#ifndef CPP_PROJECT_TEST_M_H
#define CPP_PROJECT_TEST_M_H

#include <iostream>
#include <complex>
#include <cstring>
#include <vector>

using namespace std;

template<typename T>
class matrix {
private:
    int row;

    int col;

    T *data; // data是一个包含矩阵所有元素的一维动态数组。

    int *ref_count;

public:
    int get_row();

    int get_col();

    int get_ref_count();

    int *getShape();

    matrix();

    matrix(int roe, int col);

    matrix(int row, int col, T *data); //这个constructor复数类似乎可以直接用
    matrix(const matrix<T> &other);

    matrix<T> operator=(const matrix<T> &other); // assign copy
    ~matrix();

    matrix<T> operator+(matrix<T> &other); //矩阵加法
    matrix<T> operator-(matrix<T> &other); //矩阵减法
    matrix<T> operator*(matrix<T> &other); //矩阵叉乘

    matrix<T> operator*(T number); // scalar multiplication a*2

    template<typename T1>
    friend matrix<T1> operator*(T1 number, const matrix<T1> &matrix1); // 2*a

    matrix<T> operator/(T number); // scalar division a/2

    matrix<T> transpose(); //矩阵的转置

    T dot(matrix<T> &other); //矩阵的点乘，得到一个数字
    matrix<T> reshape(int newRow, int newCol);

    //截取矩阵第二三行，第二、三列  a[1:3,1:3]
    matrix<T> slicing(int rowStart, int rowEnd, int colStart, int colEnd);

    vector<T> maxValue(string axis = "no");

    vector<T> minValue(string axis = "no");

    vector<T> sum(string axis = "no");

    vector<T> average(string axis = "no");

    matrix<T> conjugation(matrix<T> &matrix);

    void display();
};

#endif // CPP_PROJECT_TEST_M_H
