#ifndef CPP_PROJECT_TEST_M_H
#define CPP_PROJECT_TEST_M_H

#include <iostream>
#include <complex>
#include <cstring>
#include <vector>

using namespace std;

template <typename T>
class matrix
{
private:
    int row;

    int col;

    T *data; // data是一个包含矩阵所有元素的一维动态数组。

    int *ref_count;

public:
    int get_row();

    int get_col();

    T *getdata();

    int *get_ref_count();

    void setdata(int index, T a);

    T *getShape();

    matrix();

    matrix(int row, int col);

    matrix(int row, int col, T *data); //这个constructor复数类似乎可以直接用
    matrix(const matrix<T> &other);

    matrix<T> &operator=(const matrix<T> &other); // assign copy
    ~matrix();

    matrix<T> operator+(matrix<T> &other); //矩阵加法
    matrix<T> operator-(matrix<T> &other); //矩阵减法
    matrix<T> operator*(matrix<T> &other); //矩阵叉乘

    matrix<T> operator*(double number); // scalar multiplication a*2

    template <typename T1>
    friend matrix<T1> operator*(double number, matrix<T1> &matrix1); // 2*a

    matrix<T> operator/(double number); // scalar division a/2
    // v为列向量
    matrix<T> vector_multiplication(const vector<T> &v);

    // v为行向量
    template <typename T2>
    friend matrix<T2> vector_multiplication(const vector<T2> &v, const matrix<T2> &matrix1);

    matrix<T> cross(matrix<T> other);

    matrix<T> element_wise_multiplication(matrix<T> other); //矩阵对应元素相乘

    matrix<T> dot(matrix<T> &other); //矩阵的点乘,行数相同 m*1,m*n

    matrix<T> reshape(int newRow, int newCol);

    //截取矩阵第二三行，第二、三列  a[1:3,1:3]
    matrix<T> slicing(int rowStart, int rowEnd, int colStart, int colEnd);

    matrix<T> slicing(int rowStart, int colStart);

    vector<T> maxValue(string axis = "no");

    vector<T> minValue(string axis = "no");

    vector<T> sum(string axis = "no");

    vector<T> average(string axis = "no");

    matrix<T> conjugation();

    matrix<T> transpose(); //矩阵的转置

    matrix<T> convolution(matrix<T> &kernel, const string &type);

    matrix<T> rotate();

    static void caculate(T *data_temp, int temp_row, int temp_col, T *extend, const matrix<T> &kernel, int extend_col);

    void RowExchange(int a, int b);

    matrix<T> getRowEchelon() const;

    const int getRank() const;

    vector<T> EigenValues();
    vector<matrix<T>> EigenVector();
    T trace();
    matrix<T> inverse();           //实数范围求逆
    matrix<T> inverse_complex();   //复数范围求逆
    double det();                  //实数范围求行列式
    complex<double> det_complex(); //复数范围求行列式

    void display();

    void display_complexMatrix();
};

#endif // CPP_PROJECT_TEST_M_H