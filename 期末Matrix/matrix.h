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

    matrix(int row, int col);

    matrix(int row, int col, T *data); //这个constructor复数类似乎可以直接用
    matrix(const matrix<T> &other);

    matrix<T> &operator=(const matrix<T> &other); // assign copy
    ~matrix();

    matrix<T> operator+(matrix<T> &other); //矩阵加法
    matrix<T> operator-(matrix<T> &other); //矩阵减法
    matrix<T> operator*(matrix<T> &other); //矩阵叉乘

    matrix<T> operator*(T number); // scalar multiplication a*2

    template<typename T1>
    friend matrix<T1> operator*(T1 number, matrix<T1> &matrix1); // 2*a

    matrix<double> operator/(T number); // scalar division a/2
    //v为列向量
    matrix<double> vector_multiplication(const vector<T> &v);

    //v为行向量
    static matrix<T> vector_multiplication(const vector<T> &v, const matrix<T> &matrix1);

    matrix<T> cross(matrix<T> other);


    matrix<T> transpose(); //矩阵的转置

    T dot(matrix<T> &other); //矩阵的点乘，得到一个数字
    matrix<T> reshape(int newRow, int newCol);

    //截取矩阵第二三行，第二、三列  a[1:3,1:3]
    matrix<T> slicing(int rowStart, int rowEnd, int colStart, int colEnd);

    matrix<T> slicing(int rowStart, int colStart);

    vector<T> maxValue(string axis = "no");

    vector<T> minValue(string axis = "no");

    vector<T> sum(string axis = "no");

    vector<T> average(string axis = "no");

    matrix<T> conjugation(matrix<T> &matrix);//测复数

    T trace();

    double det();

    matrix<T> inverse();

    matrix<T> convolution(const matrix<T> &kernel, const string &type);

    static matrix<T> rotate(const matrix<T> &matrix1);

    void display();

    static void caculate(T *data_temp, int temp_row, int temp_col, T *extend, const matrix<T> &kernel, int extend_col);
};

#endif // CPP_PROJECT_TEST_M_H
