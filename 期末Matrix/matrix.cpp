#include "matrix.h"
#include "exception.h"
#include <iostream>
#include <complex>
#include <typeinfo>

using namespace std;

template<typename T>
void matrix<T>::display() {
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->col; ++j) {
            cout << this->data[i * this->col + j] << "\t";
        }
        cout << endl;
    }
}

template<typename T>
void display_complex(complex<T> cur) {
    if (cur.real() == 0 && cur.imag() == 0)
        cout << "0 ";
    else {
        if (cur.real())
            cout << cur.real() << "";
        if (cur.imag() > 0)
            cout << "+" << cur.imag() << "i ";
        if (cur.imag() < 0)
            cout << cur.imag() << "i ";
    }
}

template<typename T>
void matrix<T>::display_complexMatrix() {
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->col; ++j) {
            T cur = this->data[i * this->col + j];
            display_complex(cur);
            cout << "\t";
        }
        cout << endl;
    }
}

template<typename T>
int matrix<T>::get_row() {
    return this->row;
}

template<typename T>
int matrix<T>::get_col() {
    return this->col;
}

template<typename T>
T *matrix<T>::getdata() {
    return this->data;
}

template<typename T>
int *matrix<T>::get_ref_count() {
    return this->ref_count;
}

template<typename T>
T *matrix<T>::getShape() {
    return this->data;
}

template<typename T>
void matrix<T>::setdata(int index, T a) {
    this->data[index] = a;
}

template<typename T>
matrix<T>::matrix() {
//    cout << "The constructor without argument is called" << endl;
}

template<typename T>
matrix<T>::matrix(int row, int col) {
    this->row = row;
    this->col = col;
    this->data = new T[row * col]{0};
    this->ref_count = new int;
    *(this->ref_count) = 1;
//    cout << "The constructor with row and col is called, the ref_count is " << *ref_count << endl;
}

template<typename T>
matrix<T>::matrix(int row, int col, T *data) {
    this->row = row;
    this->col = col;
    this->data = new T[row * col];
    for (int i = 0; i < row * col; i++) {
        this->data[i] = data[i];
    };
    this->ref_count = new int;
    *(this->ref_count) = 1;
//    cout << "The constructor with row, col and data is called, the ref_count is " << *ref_count << endl;
}

template<typename T>
matrix<T>::matrix(const matrix<T> &other) {
    this->row = other.row;
    this->col = other.col;
    this->data = other.data;
    this->ref_count = other.ref_count;
    *(this->ref_count) += 1;
//    cout << "The copy constructor, the ref_count is " << *ref_count << endl;
}

template<typename T>
matrix<T> &matrix<T>::operator=(const matrix<T> &other) {
    this->row = other.row;
    this->col = other.col;
    *(this->ref_count) -= 1; //防止内存泄漏
    if (*(this->ref_count) == 0) {
        delete this->ref_count;
        delete[] this->data;
    }
    this->data = other.data;
    this->ref_count = other.ref_count;
    *(this->ref_count) += 1;
    return *this;
}

// avoid hard copy
template<typename T>
matrix<T>::~matrix() {
    *(this->ref_count) -= 1;
    if (*(this->ref_count) == 0) {
        delete[] this->data;
        delete this->ref_count;
    }
//    cout << "The destructor is called, the ref_count is " << *ref_count << endl;
}

//矩阵加法
template<typename T>
matrix<T> matrix<T>::operator+(matrix<T> &other) {

    if (this->row != other.row || this->col != other.col) {
        throw SizeMismatchException();
    }
    int size = (this->row) * (this->col);
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i) {
        data_temp[i] = this->data[i] + other.data[i];
    }
    return matrix(this->row, this->col, data_temp);
}

//矩阵减法
template<typename T>
matrix<T> matrix<T>::operator-(matrix<T> &other) {
    if (this->row != other.row || this->col != other.col) {
        throw SizeMismatchException();
    }
    int size = (this->row) * (this->col);
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i) {
        data_temp[i] = this->data[i] - other.data[i];
    }
    return matrix(this->row, this->col, data_temp);
}

//矩阵乘法
template<typename T>
matrix<T> matrix<T>::operator*(matrix<T> &other) {

    if (this->col != other.row) {
        throw SizeMismatchException();
    }

    int size = (this->row) * (other.col);
    T data_temp[size] = {0};
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < other.col; ++j) {
            for (int k = 0; k < this->col; ++k) {
                // cout << this->data[i * this->col + k] << " " << other.data[k * other.col + j] << endl;
                data_temp[i * other.col + j] += this->data[i * this->col + k] * other.data[k * other.col + j];
            }
        }
    }
    return matrix(this->row, other.col, data_temp);
}

// a*2 成员函数
//template<typename T>
//matrix<T> matrix<T>::operator*(double number) {
//    int size = this->row * this->col;
//    T data_temp[size] = {0};
//    for (int i = 0; i < size; ++i) {
//        data_temp[i] = this->data[i] * number;
//    }
//    return matrix(this->row, this->col, data_temp);
//}

template<typename T>
matrix<T> matrix<T>::operator*(double number) {
    int size = this->row * this->col;
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i) {
        data_temp[i] = this->data[i];
        data_temp[i] *= number;
    }
    return matrix(this->row, this->col, data_temp);
}


// 2*a
template<typename T1>
matrix<T1> operator*(double number, matrix<T1> &matrix1) {
    return matrix1 * number;
    //如何为friend函数调用constructor要像下面那样调用。
    //    return  matrix<T>(matrix1.row,matrix1.col,matrix1.temp);
}

// a/2 成员函数
template<typename T>
matrix<T> matrix<T>::operator/(double number) {
    if (number == 0) {
        throw ZeroDivideException();
    }
    int size = (this->row) * (this->col);
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i) {
        data_temp[i] = this->data[i];
        data_temp[i] /= number;
    }
    return matrix<T>(this->row, this->col, data_temp);
}

template<typename T>
matrix<T> matrix<T>::vector_multiplication(const vector<T> &v) {
    if (v.size() != this->col) {
        throw SizeMismatchException();
    }
    T data_temp[this->row] = {0};
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < v.size(); ++j) {
            data_temp[i] += this->data[i * this->col + j] * v.at(j);
        }
    }
    return matrix(this->row, 1, data_temp);
}

template<typename T2>
matrix<T2> vector_multiplication(const vector<T2> &v, const matrix<T2> &matrix1) {
    if (v.size() != matrix1.row) {
        throw SizeMismatchException();
    }
    T2 data_temp[matrix1.col];
    for (int i = 0; i < matrix1.col; ++i) {
        for (int j = 0; j < v.size(); ++j) {
            data_temp[i] += v.at(j) * matrix1.data[j * matrix1.col + i];
        }
    }
    return matrix<T2>(1, matrix1.col, data_temp);
}

template<typename T>
matrix<T> matrix<T>::cross(matrix<T> other) {
    if (this->row != 1 || this->col != 3 || other.row != 1 || other.col != 3) {
        throw Not3DVectorException();
    }
    T data_temp[3];
    data_temp[0] = this->data[1] * other.data[2] - this->data[2] * other.data[1];
    data_temp[1] = this->data[2] * other.data[0] - this->data[0] * other.data[2];
    data_temp[2] = this->data[0] * other.data[1] - this->data[1] * other.data[0];
    return matrix(1, 3, data_temp);
}

template<typename T>
matrix<T> matrix<T>::element_wise_multiplication(matrix<T> other) {
    if (this->row != other.row || this->col != other.col) {
        throw SizeMismatchException();
    }
    int size = this->row * this->col;
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i) {
        data_temp[i] = this->data[i] * other.data[i];
    }
    return matrix(this->row, this->col, data_temp);
}

template<typename T>
matrix<T> matrix<T>::dot(matrix<T> &other) {
    if (this->row != other.row || this->col != 1) {
        throw SizeMismatchException();
    }
    T data_temp[other.row * other.col] = {0};
    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < other.col; ++j) {
            data_temp[i * other.col + j] = this->data[i] * other.data[i * other.col + j];
        }
    }
    return matrix(other.row, other.col, data_temp);
}

//求矩阵转置
template<typename T>
matrix<T> matrix<T>::transpose() {
    int size = (this->row) * (this->col);
    T data_temp[size] = {0};

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            data_temp[j * this->row + i] = this->data[i * this->col + j];
        }
    }

    return matrix(this->col, this->row, data_temp);
}


template<typename T>
matrix<T> matrix<T>::reshape(int newRow, int newCol) {
    if (this->row * this->col != newRow * newCol) {
        throw InvalidReshapeSizeException();
    }
    return matrix(newRow, newCol, this->data);
}

template<typename T>
//截取矩阵第二三行，第二、三列  a[1:3,1:3]
matrix<T> matrix<T>::slicing(int rowStart, int rowEnd, int colStart, int colEnd) {
    T data_temp[(rowEnd - rowStart) * (colEnd - colStart)] = {0};
    int k = 0;
    for (int i = rowStart; i < rowEnd; ++i) {
        for (int j = colStart; j < colEnd; ++j) {
            data_temp[k++] = this->data[i * this->col + j]; // this[i][j]
        }
    }
    return matrix(rowEnd - rowStart, colEnd - colStart, data_temp);
}

template<typename T>
matrix<T> matrix<T>::slicing(int rowStart, int colStart) {
    T data_temp[(this->row - rowStart) * (this->col - colStart)] = {0};
    int k = 0;
    for (int i = rowStart; i < this->row; ++i) {
        for (int j = colStart; j < this->col; ++j) {
            data_temp[k++] = this->data[i * this->col + j]; // this[i][j]
        }
    }
    return matrix(this->row - rowStart, this->col - colStart, data_temp);
}

//求最大值
template<typename T>
vector<T> matrix<T>::maxValue(string axis) {
    if (axis != "no" && axis != "row" && axis != "col") {
        throw TypeInvalidException();
    }
    vector<T> res;
    if (axis == "no") {
        T max = this->data[0];
        for (int i = 1; i < this->row * this->col; ++i) {
            if (this->data[i] > max) {
                max = this->data[i];
            }
        }
        res.push_back(max);
    } else if (axis == "row") {
        for (int i = 0; i < this->row; i++) {
            T max = this->data[i * this->col];
            for (int j = 0; j < this->col; j++) {
                max = max < this->data[i * this->col + j] ? this->data[i * this->col + j] : max;
            }
            res.push_back(max);
        }
    } else if (axis == "col") {
        for (int i = 0; i < this->col; i++) {
            T max = this->data[col];
            for (int j = 0; j < this->row; j++) {
                max = max < this->data[j * this->col + i] ? this->data[j * this->col + i] : max;
            }
            res.push_back(max);
        }
    }
    return res;
}

//求最小值
template<typename T>
vector<T> matrix<T>::minValue(string axis) {
    if (axis != "no" && axis != "row" && axis != "col") {
        throw TypeInvalidException();
    }
    vector<T> res;
    if (axis == "no") {
        T min = this->data[0];
        for (int i = 1; i < this->row * this->col; ++i) {
            if (this->data[i] < min) {
                min = this->data[i];
            }
        }
        res.push_back(min);
    } else if (axis == "row") {
        for (int i = 0; i < this->row; i++) {
            T min = this->data[i * this->col];
            for (int j = 0; j < this->col; j++) {
                min = min > this->data[i * this->col + j] ? this->data[i * this->col + j] : min;
            }
            res.push_back(min);
        }
    } else if (axis == "col") {
        for (int i = 0; i < this->col; i++) {
            T min = this->data[col];
            for (int j = 0; j < this->row; j++) {
                min = min > this->data[j * this->col + i] ? this->data[j * this->col + i] : min;
            }
            res.push_back(min);
        }
    }
    return res;
}

//求和
template<typename T>
vector<T> matrix<T>::sum(string axis) {
    if (axis != "no" && axis != "row" && axis != "col") {
        throw TypeInvalidException();
    }
    vector<T> res;
    if (axis == "no") {
        T sum = 0;
        for (int i = 0; i < this->row * this->col; ++i) {
            sum += this->data[i];
        }
        res.push_back(sum);
    } else if (axis == "row") {
        for (int i = 0; i < this->row; i++) {
            T sum = 0;
            for (int j = 0; j < this->col; j++) {
                sum += this->data[i * this->col + j];
            }
            res.push_back(sum);
        }
    } else if (axis == "col") {
        for (int i = 0; i < this->col; i++) {
            T sum = 0;
            for (int j = 0; j < this->row; j++) {
                sum += this->data[j * this->col + i];
            }
            res.push_back(sum);
        }
    }
    return res;
}

//求平均值
template<typename T>
vector<T> matrix<T>::average(string axis) {
    if (axis != "no" && axis != "row" && axis != "col") {
        throw TypeInvalidException();
    }
    vector<T> res;
    if (axis == "no")
        res.push_back(sum()[0] /= (double) (this->row * this->col));
    else if (axis == "row") {
        res = sum("row");
        for (int i = 0; i < res.size(); i++) {
            res.at(i) /= (double) this->col;
        }
    } else if (axis == "col") {
        res = sum("col");
        for (int i = 0; i < res.size(); i++) {
            res.at(i) /= (double) this->row;
        }
    }
    return res;
}

template<typename T>
matrix<T> matrix<T>::conjugation() {
    T a;
    string x = typeid(T).name();
    if (x.substr(0, 10) == "St7complex") {
        matrix<complex<int>> res(this->row, this->col);
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                T value = this->data[i * res.col + j];
                res.data[i * res.col + j] = T(value.real(), (-1) * value.imag());
            }
        }
        return res;
    } else
        return *this;
};

//计算行列式
template<typename T>
double matrix<T>::det() {
    if (this->col != this->row)
        throw NonSquareException();

    if (this->col == 1)
        return this->data[0];
    // double *subMat = new double[(n - 1)*(n - 1)];//创建n-1阶的代数余子式阵subMat
    int n = this->col;
    matrix<T> subMat(n - 1, n - 1);
    int mov = 0;                               //判断行是否移动
    double sum = 0.0;                          // sum为行列式的值
    for (int Matrow = 0; Matrow < n; Matrow++) // Mat的行数把矩阵Mat(nn)赋值到subMat(n-1)
    {
        for (int subMatrow = 0; subMatrow < n - 1; subMatrow++) //把Mat阵第一列各元素的代数余子式存到subMat
        {
            mov = Matrow > subMatrow ? 0 : 1; // subMat中小于Matrow的行，同行赋值，等于的错过，大于的加一
            for (int j = 0; j < n - 1; j++)   //从Mat的第二列赋值到第n列
            {
                // subMat[subMatrow*(n - 1) + j] = Mat[(subMatrow + mov)*n + j + 1];
                subMat.data[subMatrow * (n - 1) + j] = this->data[(subMatrow + mov) * n + j + 1];
            }
        }
        int flag = (Matrow % 2 == 0 ? 1 : -1);
        // complex<int> flag = (Matrow % 2 == 0) ? (1,0) : (-1,0); //因为列数为0，所以行数是偶数时候，代数余子式为1.
        // sum += flag* Mat[Matrow*n] * det(n - 1, subMat);//Mat第一列各元素与其代数余子式积的和即为行列式
        sum += flag * this->data[Matrow * n] * subMat.det();
    }
    // delete subMat;
    return sum;
}

template<typename T>
complex<double> matrix<T>::det_complex() {
    if (this->col != this->row)
        throw NonSquareException();

    if (this->col == 1)
        return this->data[0];
    // double *subMat = new double[(n - 1)*(n - 1)];//创建n-1阶的代数余子式阵subMat
    int n = this->col;
    matrix<T> subMat(n - 1, n - 1);
    int mov = 0;                               //判断行是否移动
    complex<double> sum(0, 0);                          // sum为行列式的值
    for (int Matrow = 0; Matrow < n; Matrow++) // Mat的行数把矩阵Mat(nn)赋值到subMat(n-1)
    {
        for (int subMatrow = 0; subMatrow < n - 1; subMatrow++) //把Mat阵第一列各元素的代数余子式存到subMat
        {
            mov = Matrow > subMatrow ? 0 : 1; // subMat中小于Matrow的行，同行赋值，等于的错过，大于的加一
            for (int j = 0; j < n - 1; j++)   //从Mat的第二列赋值到第n列
            {
                // subMat[subMatrow*(n - 1) + j] = Mat[(subMatrow + mov)*n + j + 1];
                subMat.data[subMatrow * (n - 1) + j] = this->data[(subMatrow + mov) * n + j + 1];
            }
        }
        complex<double> positive(1, 0);
        complex<double> negetive(-1, 0);
        complex<double> flag;
        if (Matrow % 2 == 0) flag = positive;
        else flag = negetive;
//        complex<double> flag = () ? (1, 0) : (-1, 0);
        // sum += flag* Mat[Matrow*n] * det(n - 1, subMat);//Mat第一列各元素与其代数余子式积的和即为行列式
        sum += flag * this->data[Matrow * n] * subMat.det_complex();
    }
    // delete subMat;
    return sum;
}


template<typename T>
matrix<T> matrix<T>::inverse() {
    int n = this->col;
    matrix<T> subMat(n - 1, n - 1);
    matrix<T> adjointMat(n, n); // Mat的伴随矩阵adjointMat[]
    T detMat = this->det();     // Mat的行列式
    matrix<T> inverseMat(n, n); // Mat的逆矩阵inverseMat[]

    if (detMat == 0)
        throw IrreversibleException();
    else {
        int adjointMat_index = 0;
        for (int Mat_i = 0; Mat_i < n; Mat_i++) {
            for (int Mat_j = 0; Mat_j < n; Mat_j++) {
                int M_index = 0;
                //求第Mat_i,Mat_j的余子式Mij
                for (int i = 0; i < n * n; i++) //循环Mat的下标选出余子式
                {
                    //第x个元素对应n阶矩阵的第row行第col列
                    int x = i;
                    int row = 0;
                    int col = 0;
                    while (1) {
                        if (x - n >= 0) {
                            x = x - n;
                            row++;
                        } else {
                            col = x;
                            break;
                        }
                    }
                    if (row != Mat_i && col != Mat_j) {
                        subMat.data[M_index] = this->data[i];
                        M_index++;
                    }
                }
                adjointMat.data[adjointMat_index] = ((Mat_i + Mat_j) % 2 == 0 ? 1 : -1) * subMat.det(); //求伴随矩阵各元素
                adjointMat_index++;
            }
        }
        // adjointMat = MatT(n, n, adjointMat); //转置
        adjointMat = adjointMat.transpose();

        cout << "---- The adjoin matrix is ----" << endl;
        int adMat_index = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << adjointMat.data[adMat_index] << " ";
                adMat_index++;
            }
            cout << endl;
        }

        for (int i = 0; i < n * n; i++) {
            inverseMat.setdata(i, adjointMat.data[i] / detMat);
        }
        cout << "---- The reverse matrix is ----" << endl;
        int inverseMat_index = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << inverseMat.data[inverseMat_index] << "\t";
                inverseMat_index++;
            }
            cout << endl;
        }
    }
    return inverseMat;
}

template<typename T>
matrix<T> matrix<T>::inverse_complex() {
    int n = this->col;
    matrix<T> subMat(n - 1, n - 1);
    matrix<T> adjointMat(n, n); // Mat的伴随矩阵adjointMat[]
    T detMat = this->det_complex();     // Mat的行列式
    matrix<T> inverseMat(n, n); // Mat的逆矩阵inverseMat[]

    complex<double> zero(0, 0);
    if (detMat == zero)
        throw IrreversibleException();
    else {
        int adjointMat_index = 0;
        for (int Mat_i = 0; Mat_i < n; Mat_i++) {
            for (int Mat_j = 0; Mat_j < n; Mat_j++) {
                int M_index = 0;
                //求第Mat_i,Mat_j的余子式Mij
                for (int i = 0; i < n * n; i++) //循环Mat的下标选出余子式
                {
                    //第x个元素对应n阶矩阵的第row行第col列
                    int x = i;
                    int row = 0;
                    int col = 0;
                    while (1) {
                        if (x - n >= 0) {
                            x = x - n;
                            row++;
                        } else {
                            col = x;
                            break;
                        }
                    }
                    if (row != Mat_i && col != Mat_j) {
                        subMat.data[M_index] = this->data[i];
                        M_index++;
                    }
                }
                complex<double> positive(1, 0);
                complex<double> negetive(-1, 0);
                complex<double> flag;
                if ((Mat_i + Mat_j) % 2 == 0) flag = positive;
                else flag = negetive;
                adjointMat.data[adjointMat_index] = flag * subMat.det_complex(); //求伴随矩阵各元素
                adjointMat_index++;
            }
        }
        // adjointMat = MatT(n, n, adjointMat); //转置
        adjointMat = adjointMat.transpose();

        cout << "---- The adjoin matrix is ----" << endl;
        int adMat_index = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << adjointMat.data[adMat_index] << " ";
                adMat_index++;
            }
            cout << endl;
        }

        for (int i = 0; i < n * n; i++) {
            inverseMat.setdata(i, adjointMat.data[i] / detMat);
        }
        cout << "---- The reverse matrix is ----" << endl;
        int inverseMat_index = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << inverseMat.data[inverseMat_index] << "\t";
                inverseMat_index++;
            }
            cout << endl;
        }
    }
    return inverseMat;
}

template<typename T>
T matrix<T>::trace() {
    if (this->row != this->col) {
        throw NonSquareException();
    }
    T cnt = 0;
    for (int i = 0; i < this->row; ++i) {
        cnt += this->data[i * this->col + i];
    }
    return cnt;
}

//把矩阵旋转180度
template<typename T>
matrix<T> matrix<T>::rotate() {
    T data_temp[this->row * this->col];
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->col; j++) {
            data_temp[i * this->col + j] = this->data[(this->row - 1 - i) * this->col + this->col - 1 - j];
        }
    }
    return matrix<T>(this->row, this->col, data_temp);
}

//卷积运算
template<typename T>
matrix<T> matrix<T>::convolution(matrix<T> &in, const string &type) {
    if (type != "full" && type != "valid" && type != "same") {
        throw TypeInvalidException();
    }
    matrix<T> kernel = in.rotate();
    int temp_row;
    int temp_col;
    int size;
    if (type == "valid") {
        // valid 型 不扩充
        temp_row = this->row - kernel.row + 1;
        temp_col = this->col - kernel.row + 1;
        size = (temp_row) * (temp_col);
        T data_temp[size] = {0};
        caculate(data_temp, temp_row, temp_col, this->data, kernel, this->col);
        return matrix(temp_row, temp_col, data_temp);
    } else if (type == "same") {
        // same 型 扩充一圈
        T extend[(this->row + 2) * (this->col + 2)] = {0}; //扩充初始矩阵，填0
        for (int i = 1; i < this->row + 1; ++i) {
            for (int j = 1; j < this->col + 1; ++j) {
                extend[i * (this->col + 2) + j] = this->data[(i - 1) * (this->col) + (j - 1)];
            }
        }
        temp_row = this->row - kernel.row + 3;
        temp_col = this->col - kernel.row + 3;
        size = (temp_row) * (temp_col);
        T data_temp[size] = {0};
        caculate(data_temp, temp_row, temp_col, extend, kernel, this->col + 2);
        return matrix(temp_row, temp_col, data_temp);
    } else if (type == "full") {
        //扩充两圈
        T extend[(this->row + 4) * (this->col + 4)] = {0}; //扩充初始矩阵，填0
        for (int i = 2; i < this->row + 2; ++i) {
            for (int j = 2; j < this->col + 2; ++j) {
                extend[i * (this->col + 4) + j] = this->data[(i - 2) * (this->col) + (j - 2)];
            }
        }
        temp_row = this->row - kernel.row + 5;
        temp_col = this->col - kernel.row + 5;
        size = (temp_row) * (temp_col);
        T data_temp[size] = {0};
        caculate(data_temp, temp_row, temp_col, extend, kernel, this->col + 4);
        return matrix(temp_row, temp_col, data_temp);
    }
}

//卷积操作
template<typename T>
void
matrix<T>::caculate(T *data_temp, int temp_row, int temp_col, T *extend, const matrix<T> &kernel, int extend_col) {
    for (int i = 0; i < temp_row; ++i) {
        for (int j = 0; j < temp_col; ++j) {
            T sum = 0;
            for (int m = 0; m < kernel.row; ++m) {
                for (int n = 0; n < kernel.col; ++n) {
                    sum = sum + extend[(i + m) * extend_col + j + n] * kernel.data[m * kernel.col + n];
                }
            }
            data_temp[i * temp_col + j] = sum;
        }
    }
}

template<typename T>
vector<T> matrix<T>::EigenValues() {
    vector<T> res;
    int MIN = -100;
    int MAX = 100;
    // int sign = 0;
    double STEP = 1;
    double PRE = 0.0001;

    /*复制一份矩阵参与运算。*/
    matrix<T> c_mat_alt(this->row, this->col, this->data);
    double real;
    for (real = MIN; real <= MAX; real += STEP) {
        for (int i = 0; i < this->row * this->col; i++) {
            c_mat_alt.data[i] = this->data[i];
        }
        for (int i = 0; i < this->col; i++) {
            c_mat_alt.data[i * this->col + i] -= real;
        }
        /*如果行列式计算结果显示实部和虚部均小于给定精度（被认作为0），判定real + imag i是特征值。*/
        if (fabs(c_mat_alt.det()) <= PRE) {
            res.push_back(real);
        }
    }
    return res;
}

template<typename T>
vector<matrix<T>> matrix<T>::EigenVector() {
    vector<T> eigenValue = this->EigenValues();
    vector<matrix<T>> values;
    matrix<T> c_ALU(this->row, this->col, this->data);

    for (int j = 0; j < eigenValue.size(); j++) {
        T cur = eigenValue.at(j);
        for (int i = 0; i < this->row * this->col; i++) {
            c_ALU.data[i] = this->data[i];
        }
        for (int i = 0; i < this->col; i++) {
            c_ALU.data[i * this->col + i] -= cur;
        }
        c_ALU = c_ALU.getRowEchelon();
        int rank = c_ALU.getRank();
        int numberOfValue = this->row - rank;
        int v = 0; //有n列无主元

        for (int jj = 0; jj < c_ALU.col; jj++) {
            for (int ii = jj - v; ii < c_ALU.row; ii++) { //对于第[j]列 本来是第[j]行 但因为跳过了v 我们的主元要从第[j-v]行开始找主元
                if (c_ALU.data[ii * this->col + jj] == 0) { //当前行主元列为0
                    if (ii == row - 1) { //整列都是0 没有主元
                        T arry[c_ALU.col]{0};
                        for (int l = 0; l < c_ALU.row; l++) {
                            if (c_ALU.data[l * c_ALU.col + jj] != 0)
                                arry[l] = (-1) * c_ALU.data[l * c_ALU.col + jj];
                        }
                        arry[jj] = 1;
                        values.push_back(matrix(this->col, 1, arry));
                        v++;
                        break;
                    } else {
                        continue;
                    }
                } else
                    break;
            }
        }
    }
    return values;
}


template<typename T>
const int matrix<T>::getRank() const {
    matrix<T> Matrix_RE = this->getRowEchelon();
    //得到简化阶梯矩阵来求秩
    // Matrix_RE.show();

    int rows = Matrix_RE.row;
    int cols = Matrix_RE.col;
    int rank = rows;
    for (int i = 0; i < rows; i++) {
        bool allzero = true;
        for (int j = 0; j < cols; j++) {
            if (Matrix_RE.data[i * this->col + j] != 0) {
                allzero = false;
                break;
            }
        }
        if (allzero) {
            rank--;
        }
    }
    return rank;
}

template<typename T>
matrix<T> matrix<T>::getRowEchelon() const {
    matrix<T> Matrix_RE = matrix<T>(*this);
    int v = 0; //有n列无主元
    int factor = 1;
    for (int j = 0; j < Matrix_RE.col; j++) {
        for (int i = j - v; i < Matrix_RE.row; i++) { //对于第[j]列 本来是第[j]行 但因为跳过了v 我们的主元要从第[j-v]行开始找主元
            if (Matrix_RE.data[i * this->col + j] == 0) { //当前行主元列为0
                if (i == row - 1) { //整列都是0 没有主元
                    v++;
                    break;
                } else {
                    continue;
                }
            } else { //找到非零元 先交换到第j-v行 然后把其他都化为0
                Matrix_RE.RowExchange(i, j - v);
                // T pivot = Matrix_RE.vectors[j - v][j];//找到首元位置
                T pivot = Matrix_RE.data[(j - v) * this->col + j];
                for (int i_cur = 0; i_cur < row; i_cur++) {
                    T AtPivotCol = Matrix_RE.data[(i_cur) * this->col + j];
                    // T AtPivotCol = Matrix_RE.vectors[i_cur][j];//看第i_cur行的主元列是不是0 如果不是 每一项都乘以主元然后和主元*AtPivotCol相减成0
                    if (AtPivotCol != 0 && i_cur != j - v) {
                        for (int j_cur = 0; j_cur < col; j_cur++) {
                            Matrix_RE.data[(i_cur) * this->col + j_cur] *= pivot;
                            Matrix_RE.data[(j - v) * this->col + j_cur] *= AtPivotCol;
                            Matrix_RE.data[(i_cur) * this->col + j_cur] -= Matrix_RE.data[(j - v) * this->col +
                                                                                          j_cur];
                            // Matrix_RE.vectors[i_cur][j_cur] /= pivot;
                            Matrix_RE.data[(j - v) * this->col + j_cur] /= AtPivotCol;
                        }
                    }
                }
                for (int o = 0; o < Matrix_RE.col; o++) {
                    Matrix_RE.data[(j - v) * this->col + o] /= pivot;
                }
                break;
            }
        }
    }

    return Matrix_RE;
}

template<typename T>
void matrix<T>::RowExchange(int a, int b) {
    T temp[this->col];
    for (int i = 0; i < this->col; i++) {
        temp[i] = this->data[a * this->col + i];
    }
    for (int i = 0; i < this->col; i++) {
        this->data[a * this->col + i] = this->data[b * this->col + i];
    }
    for (int i = 0; i < this->col; i++) {
        this->data[b * this->col + i] = temp[i];
    }
}