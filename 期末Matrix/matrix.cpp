#include "matrix.h"
#include "exception.h"
#include <iostream>
#include <complex>

using namespace std;

template<typename T>
void matrix<T>::display() {
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->col; ++j) {
            cout << this->data[i * this->col + j] << " ";
        }
        cout << endl;
    }
}

template<typename T>
int matrix<T>::get_row() {
    return this->nrows;
}

template<typename T>
int matrix<T>::get_col() {
    return this->get_col();
}

template<typename T>
int matrix<T>::get_ref_count() {
    return this->ref_count;
}

template<typename T>
int *matrix<T>::getShape() {
    return this->data;
}

template<typename T>
matrix<T>::matrix() {}

template<typename T>
matrix<T>::matrix(int row, int col) {
    this->row = row;
    this->col = col;
    this->data = new T[row * col]{0};
    this->ref_count = new int;
    *(this->ref_count) = 1;
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
}

template<typename T>
matrix<T>::matrix(const matrix<T> &other) {
    this->row = other.row;
    this->col = other.col;
    this->data = other.data;
    this->ref_count = other.ref_count;
    *(this->ref_count) += 1;
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
template<typename T>
matrix<T> matrix<T>::operator*(T number) {
    int size = this->row * this->col;
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i) {
        data_temp[i] = this->data[i] * number;
    }
    return matrix(this->row, this->col, data_temp);
}

// 2*a
template<typename T1>
matrix<T1> operator*(T1 number, matrix<T1> &matrix1) {
    return matrix1 * number;
    //如何为friend函数调用constructor要像下面那样调用。
    //    return  matrix<T>(matrix1.row,matrix1.col,matrix1.temp);
}

// a/2 成员函数
template<typename T>
matrix<double> matrix<T>::operator/(T number) {
    if (number == 0) {
        throw ZeroDivideException();
    }
    int size = this->row * this->col;
    double data_temp[size] = {0};
    for (int i = 0; i < size; ++i) {
        data_temp[i] = (double) this->data[i] / number;
    }
    return matrix<double>(this->row, this->col, data_temp);
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
T matrix<T>::dot(matrix<T> &other) {
    if (this->row != other.row || this->col != other.col) {
        throw SizeMismatchException();
    }
    int size = this->row * this->col;
    T dot_answer;
    for (int i = 0; i < size; ++i) {
        dot_answer += this->data[i] * other.data[i];
    }
    return dot_answer;
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
    vector<T> res;
    if (axis == "no")
        res.push_back((double) sum()[0] / (this->row * this->col));
    else if (axis == "row") {
        res = sum("row");
        for (int i = 0; i < res.size(); i++) {
            res.at(i) = (double) res.at(i) / this->col;
        }
    } else if (axis == "col") {
        res = sum("col");
        for (int i = 0; i < res.size(); i++) {
            res.at(i) = (double) res.at(i) / this->row;
        }
    }
    return res;
}

template<typename T>
matrix<T> matrix<T>::conjugation(matrix<T> &original) {
    T a;
    string x = typeid(T).name();
    if (x.substr(0, 10) == "St7complex") {
        matrix<T> res(original.row, original.col);
        for (int i = 0; i < original.row; i++) {
            for (int j = 0; j < original.col; j++) {
                complex<T> value = original.data[i * res.col + j];
                res.data[i * res.col + j] = complex<T>(value.real, (-1) * value.imag);
            }
        }
        return res;
    } else
        return original;
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
        int flag = (Matrow % 2 == 0 ? 1 : -1); //因为列数为0，所以行数是偶数时候，代数余子式为1.
        // sum += flag* Mat[Matrow*n] * det(n - 1, subMat);//Mat第一列各元素与其代数余子式积的和即为行列式
        sum += flag * this->data[Matrow * n] * subMat.det();
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

    // double *subMat = new double[(n - 1)*(n - 1)];//Mat的(代数)余子式subMat[]
    // double *adjointMat = new double[n*n];
    // double detMat;
    // double *inverseMat = new double[n*n];
    // detMat = det(n, Mat);//求Mat的行列式
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

//        cout << "该矩阵的伴随矩阵为" << endl;
//        int adMat_index = 0;
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                cout << adjointMat.data[adMat_index] << "\t";
//                adMat_index++;
//            }
//            cout << endl;
//        }

        for (int i = 0; i < n * n; i++)
            inverseMat.data[i] = adjointMat.data[i] / detMat;
        // cout << "该矩阵的逆矩阵为" << endl;
        // int inverseMat_index = 0;
        // for (i = 0; i < n; i++)
        // {
        //     for (j = 0; j < n; j++)
        //     {
        //         cout << inverseMat[inverseMat_index] << "\t";
        //         inverseMat_index++;
        //     }
        //     cout << endl;
        // }
    }
    return inverseMat;
}#include "matrix.h"
#include "exception.h"
#include <iostream>
#include <complex>

using namespace std;

template<typename T>
void matrix<T>::display() {
    for (int i = 0; i < this->row; i++) {
        for (int j = 0; j < this->col; ++j) {
            cout << this->data[i * this->col + j] << " ";
        }
        cout << endl;
    }
}

template<typename T>
int matrix<T>::get_row() {
    return this->nrows;
}

template<typename T>
int matrix<T>::get_col() {
    return this->get_col();
}

template<typename T>
int matrix<T>::get_ref_count() {
    return this->ref_count;
}

template<typename T>
int *matrix<T>::getShape() {
    return this->data;
}

template<typename T>
matrix<T>::matrix() {}

template<typename T>
matrix<T>::matrix(int row, int col) {
    this->row = row;
    this->col = col;
    this->data = new T[row * col]{0};
    this->ref_count = new int;
    *(this->ref_count) = 1;
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
}

template<typename T>
matrix<T>::matrix(const matrix<T> &other) {
    this->row = other.row;
    this->col = other.col;
    this->data = other.data;
    this->ref_count = other.ref_count;
    *(this->ref_count) += 1;
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
template<typename T>
matrix<T> matrix<T>::operator*(T number) {
    int size = this->row * this->col;
    T data_temp[size] = {0};
    for (int i = 0; i < size; ++i) {
        data_temp[i] = this->data[i] * number;
    }
    return matrix(this->row, this->col, data_temp);
}

// 2*a
template<typename T1>
matrix<T1> operator*(T1 number, matrix<T1> &matrix1) {
    return matrix1 * number;
    //如何为friend函数调用constructor要像下面那样调用。
    //    return  matrix<T>(matrix1.row,matrix1.col,matrix1.temp);
}

// a/2 成员函数
template<typename T>
matrix<double> matrix<T>::operator/(T number) {
    if (number == 0) {
        throw ZeroDivideException();
    }
    int size = this->row * this->col;
    double data_temp[size] = {0};
    for (int i = 0; i < size; ++i) {
        data_temp[i] = (double) this->data[i] / number;
    }
    return matrix<double>(this->row, this->col, data_temp);
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
T matrix<T>::dot(matrix<T> &other) {
    if (this->row != other.row || this->col != other.col) {
        throw SizeMismatchException();
    }
    int size = this->row * this->col;
    T dot_answer;
    for (int i = 0; i < size; ++i) {
        dot_answer += this->data[i] * other.data[i];
    }
    return dot_answer;
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
    vector<T> res;
    if (axis == "no")
        res.push_back((double) sum()[0] / (this->row * this->col));
    else if (axis == "row") {
        res = sum("row");
        for (int i = 0; i < res.size(); i++) {
            res.at(i) = (double) res.at(i) / this->col;
        }
    } else if (axis == "col") {
        res = sum("col");
        for (int i = 0; i < res.size(); i++) {
            res.at(i) = (double) res.at(i) / this->row;
        }
    }
    return res;
}

template<typename T>
matrix<T> matrix<T>::conjugation(matrix<T> &original) {
    T a;
    string x = typeid(T).name();
    if (x.substr(0, 10) == "St7complex") {
        matrix<T> res(original.row, original.col);
        for (int i = 0; i < original.row; i++) {
            for (int j = 0; j < original.col; j++) {
                complex<T> value = original.data[i * res.col + j];
                res.data[i * res.col + j] = complex<T>(value.real, (-1) * value.imag);
            }
        }
        return res;
    } else
        return original;
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
        int flag = (Matrow % 2 == 0 ? 1 : -1); //因为列数为0，所以行数是偶数时候，代数余子式为1.
        // sum += flag* Mat[Matrow*n] * det(n - 1, subMat);//Mat第一列各元素与其代数余子式积的和即为行列式
        sum += flag * this->data[Matrow * n] * subMat.det();
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

    // double *subMat = new double[(n - 1)*(n - 1)];//Mat的(代数)余子式subMat[]
    // double *adjointMat = new double[n*n];
    // double detMat;
    // double *inverseMat = new double[n*n];
    // detMat = det(n, Mat);//求Mat的行列式
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

//        cout << "该矩阵的伴随矩阵为" << endl;
//        int adMat_index = 0;
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                cout << adjointMat.data[adMat_index] << "\t";
//                adMat_index++;
//            }
//            cout << endl;
//        }

        for (int i = 0; i < n * n; i++)
            inverseMat.data[i] = adjointMat.data[i] / detMat;
        // cout << "该矩阵的逆矩阵为" << endl;
        // int inverseMat_index = 0;
        // for (i = 0; i < n; i++)
        // {
        //     for (j = 0; j < n; j++)
        //     {
        //         cout << inverseMat[inverseMat_index] << "\t";
        //         inverseMat_index++;
        //     }
        //     cout << endl;
        // }
    }
    return inverseMat;
}
