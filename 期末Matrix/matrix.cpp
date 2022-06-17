#include "matrix.h"
#include "exception.h"
#include <iostream>
#include <complex>
using namespace std;

template <typename T>
void matrix<T>::display()
{
    for (int i = 0; i < this->row; i++)
    {
        for (int j = 0; j < this->col; ++j)
        {
            cout << this->data[i * this->col + j] << " ";
        }
        cout << endl;
    }
}

template <typename T>
int matrix<T>::get_row()
{
    return this->nrows;
}

template <typename T>
int matrix<T>::get_col()
{
    return this->get_col();
}

template <typename T>
int matrix<T>::get_ref_count()
{
    return this->ref_count;
}

template <typename T>
int *matrix<T>::getShape()
{
    return this->data;
}

template <typename T>
matrix<T>::matrix() {}

template <typename T>
matrix<T>::matrix(int row, int col)
{
    this->row = row;
    this->col = col;
    this->data = new T[row * col]{0};
    *(this->ref_count) = 1;
}

template <typename T>
matrix<T>::matrix(int row, int col, T *data)
{
    this->row = row;
    this->col = col;
    this->data = new T[row * col];
    for (int i = 0; i < row * col; i++)
    {
        this->data[i] = data[i];
    };
    this->ref_count = new int;
    *(this->ref_count) = 1;
}

template <typename T>
matrix<T>::matrix(const matrix<T> &other)
{
    this->row = other.row;
    this->col = other.col;
    this->data = other.data;
    this->ref_count = other.ref_count;
    *(this->ref_count) += 1;
}

template <typename T>
matrix<T> matrix<T>::operator=(const matrix<T> &other)
{
    this->row = other.row;
    this->col = other.col;
    *(this->ref_count) -= 1; //防止内存泄漏
    if (*(this->ref_count) == 0 && this->data == nullptr)
    {
        delete this->ref_count;
        delete[] this->data;
    }
    this->data = other.data;
    this->ref_count = other.ref_count;
    *(this->ref_count) += 1;
    return *this;
}

// avoid hard copy
template <typename T>
matrix<T>::~matrix()
{
    *(this->ref_count) -= 1;
    if (*(this->ref_count) == 0 && this->data == nullptr)
    {
        delete[] this->data;
        delete this->ref_count;
    }
}

//释放内存？？？？？？？？？？？
//矩阵加法
template <typename T>
matrix<T> matrix<T>::operator+(matrix<T> &other)
{

    if (this->row != other.row || this->col != other.col)
    {
        throw SizeMismatchException();
    }
    int size = (this->row) * (this->col);
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i)
    {
        data_temp[i] = this->data[i] + other.data[i];
    }
    return matrix(this->row, this->col, data_temp);
}

//矩阵减法
template <typename T>
matrix<T> matrix<T>::operator-(matrix<T> &other)
{
    if (this->row != other.row || this->col != other.col)
    {
        throw SizeMismatchException();
    }
    int size = (this->row) * (this->col);
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i)
    {
        data_temp[i] = this->data[i] - other.data[i];
    }
    return matrix(this->row, this->col, data_temp);
}

//矩阵乘法
template <typename T>
matrix<T> matrix<T>::operator*(matrix<T> &other)
{

    if (this->col != other.row)
    {
        throw SizeMismatchException();
    }

    int size = (this->row) * (other.col);
    T data_temp[size] = {0};
    for (int i = 0; i < this->row; ++i)
    {
        for (int j = 0; j < other.col; ++j)
        {
            for (int k = 0; k < this->col; ++k)
            {
                cout << this->data[i * this->col + k] << " " << other.data[k * other.col + j] << endl;
                data_temp[i * other.col + j] += this->data[i * this->col + k] * other.data[k * other.col + j];
            }
        }
    }
    return matrix(this->row, other.col, data_temp);
}

// a*2 成员函数
template <typename T>
matrix<T> matrix<T>::operator*(T number)
{
    int size = this->row * this->col;
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i)
    {
        data_temp[i] = this->data[i] * number;
    }
    return matrix(this->row, this->col, data_temp);
}

// 2*a
template <typename T>
matrix<T> operator*(T number, matrix<T> &matrix1)
{
    return matrix1 * number;
    //如何为friend函数调用constructor要像下面那样调用。
    //    return  matrix<T>(matrix1.row,matrix1.col,matrix1.temp);
}

// a/2 成员函数
template <typename T>
matrix<T> matrix<T>::operator/(T number)
{
    if (number == 0)
    {
        throw ZeroDivideException();
    }
    int size = this->row * this->col;
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i)
    {
        data_temp[i] = this->data[i] / number;
    }
    return matrix(this->row, this->col, data_temp);
}

//求矩阵转置
template <typename T>
matrix<T> matrix<T>::transpose()
{
    int size = (this->row) * (this->col);
    T data_temp[size] = {0};
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            data_temp[j * this->row + i] = this->data[i * this->col + j];
        }
    }
}

template <typename T>
T matrix<T>::dot(matrix<T> &other)
{
    if (this->row != other.row || this->col != other.col)
    {
        throw SizeMismatchException();
    }
    int size = this->row * this->col;
    T dot_answer;
    for (int i = 0; i < size; ++i)
    {
        dot_answer += this->data[i] * other.data[i];
    }
    return dot_answer;
}

template <typename T>
matrix<T> matrix<T>::reshape(int newRow, int newCol)
{
    if (this->row * this->col != newRow * newCol)
    {
        throw InvalidReshapeSizeException();
    }
    matrix<T> reMatrix(newRow, newCol, this->data); //??结果生成了又不见了
}

template <typename T> //截取矩阵第二三行，第二、三列  a[1:3,1:3]
matrix<T> matrix<T>::slicing(int rowStart, int rowEnd, int colStart, int colEnd)
{
    T data_temp[(rowEnd - rowStart) * (colEnd - colStart)] = {0};
    int k = 0;
    for (int i = rowStart; i < rowEnd; ++i)
    {
        for (int j = colStart; j < colEnd; ++j)
        {
            data_temp[k++] = this->data[i * this->col + j]; // this[i][j]
        }
    }
    return matrix(rowEnd - rowStart, colEnd - colStart, data_temp);
}

//求最大值
template <typename T>
vector<T> matrix<T>::maxValue(string axis)
{
    vector<T> res;
    if (axis == "no")
    {
        T max = this->data[0];
        for (int i = 1; i < this->row * this->col; ++i)
        {
            if (this->data[i] > max)
            {
                max = this->data[i];
            }
        }
        res.push_back(max);
    }
    else if (axis == "row")
    {
        for (int i = 0; i < this->row; i++)
        {
            T max = this->data[i * this->col];
            for (int j = 0; j < this->col; j++)
            {
                max = max < this->data[i * this->col + j] ? this->data[i * this->col + j] : max;
            }
            res.push_back(max);
        }
    }
    else if (axis == "col")
    {
        for (int i = 0; i < this->col; i++)
        {
            T max = this->data[col];
            for (int j = 0; j < this->row; j++)
            {
                max = max < this->data[j * this->col + i] ? this->data[j * this->col + i] : max;
            }
            res.push_back(max);
        }
    }
    return max;
}

//求最小值
template <typename T>
vector<T> matrix<T>::minValue(string axis)
{
    vector<T> res;
    if (axis == "no")
    {
        T min = this->data[0];
        for (int i = 1; i < this->row * this->col; ++i)
        {
            if (this->data[i] < min)
            {
                min = this->data[i];
            }
        }
        res.push_back(min);
    }
    else if (axis == "row")
    {
        for (int i = 0; i < this->row; i++)
        {
            T min = this->data[i * this->col];
            for (int j = 0; j < this->col; j++)
            {
                min = min > this->data[i * this->col + j] ? this->data[i * this->col + j] : min;
            }
            res.push_back(min);
        }
    }
    else if (axis == "col")
    {
        for (int i = 0; i < this->col; i++)
        {
            T min = this->data[col];
            for (int j = 0; j < this->row; j++)
            {
                min = min > this->data[j * this->col + i] ? this->data[j * this->col + i] : min;
            }
            res.push_back(min);
        }
    }
    return min;
}

//求和
template <typename T>
vector<T> matrix<T>::sum(string axis)
{
    vector<T> res;
    if (axis = "no")
    {
        T sum = 0;
        for (int i = 0; i < this->row * this->col; ++i)
        {
            sum += this->data[i];
        }
        res.push_back(sum);
    }
    else if (axis == "row")
    {
        for (int i = 0; i < this->row; i++)
        {
            T sum = 0;
            for (int j = 0; j < this->col; j++)
            {
                sum += this->data[i * this->col + j];
            }
            res.push_back(sum);
        }
    }
    else if (axis == "col")
    {
        for (int i = 0; i < this->col; i++)
        {
            T sum = 0;
            for (int j = 0; j < this->row; j++)
            {
                sum += this->data[j * this->col + i];
            }
            res.push_back(sum);
        }
    }
    return res;
}

//求平均值
template <typename T>
vector<T> matrix<T>::average(string axis)
{
    vector<T> res;
    if (axis = "no")
        res.push_back((double)sum() / (this->row * this->col));
    else if (axis == "row")
    {
        res = sum("row");
        for (int i = 0; i < res.size(); i++)
        {
            res.at(i) = (double)res.at(i) / this->col;
        }
    }
    else if (axis == "col")
    {
        res = sum("col");
        for (int i = 0; i < res.size(); i++)
        {
            res.at(i) = (double)res.at(i) / this->row;
        }
    }
    return res;
}

template <typename T>
matrix<T> matrix<T>::conjugation(matrix<T> &original)
{
    T a;
    string x = typeid(T).name();
    if (x.substr(0, 10) == "St7complex")
    {
        matrix<T> res(original.row, original.col);
        for (int i = 0; i < original.row; i++)
        {
            for (int j = 0; j < original.col; j++)
            {
                complex<T> value = original.data[i * res.col + j];
                res.data[i * res.col + j] = complex<T>(value.real, (-1) * value.imag);
            }
        }
        return res;
    }
    else
        return original;
};
