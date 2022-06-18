#include "matrix.h"
#include "exception.h"
#include <iostream>
#include <complex>
#include <math.h>
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
    this->ref_count = new int;
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
    if (*(this->ref_count) == 0)
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
    if (*(this->ref_count) == 0)
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
                // cout << this->data[i * this->col + k] << " " << other.data[k * other.col + j] << endl;
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
template <typename T1>
matrix<T1> operator*(T1 number, matrix<T1> &matrix1)
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

    return matrix(this->col, this->row, data_temp);
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

//计算行列式
template <typename T>
double matrix<T>::det()
{
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

template <typename T>
matrix<T> matrix<T>::inverse()
{
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
        cout << "该矩阵不可逆!" << endl;
    else
    {
        int adjointMat_index = 0;
        for (int Mat_i = 0; Mat_i < n; Mat_i++)
        {
            for (int Mat_j = 0; Mat_j < n; Mat_j++)
            {
                int M_index = 0;
                //求第Mat_i,Mat_j的余子式Mij
                for (int i = 0; i < n * n; i++) //循环Mat的下标选出余子式
                {
                    //第x个元素对应n阶矩阵的第row行第col列
                    int x = i;
                    int row = 0;
                    int col = 0;
                    while (1)
                    {
                        if (x - n >= 0)
                        {
                            x = x - n;
                            row++;
                        }
                        else
                        {
                            col = x;
                            break;
                        }
                    }
                    if (row != Mat_i && col != Mat_j)
                    {
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

        cout << "该矩阵的伴随矩阵为" << endl;
        int adMat_index = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cout << adjointMat.data[adMat_index] << "\t";
                adMat_index++;
            }
            cout << endl;
        }

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

template <typename T>
vector<T> matrix<T>::EigenValues()
{
    vector<T> res;
    int MIN = 0;
    int MAX = 10;
    // int sign = 0;
    double STEP = 0.001;
    double PRE = 0.0001;

    /*复制一份矩阵参与运算。*/
    matrix<T> c_mat_alt(this->row, this->col, this->data);
    double real;
    for (real = MIN; real <= MAX; real += STEP)
    {
        for (int i = 0; i < this->row + this->col; i++)
        {
            c_mat_alt.data[i] = this->data[i];
        }
        for (int i = 0; i < this->col; i++)
        {
            c_mat_alt.data[i * this->col + i] -= real;
        }
        /*如果行列式计算结果显示实部和虚部均小于给定精度（被认作为0），判定real + imag i是特征值。*/
        if (fabs(c_mat_alt.det()) <= PRE)
        {
            res.push_back(real);
        }
    }
    return res;
}

// template <typename T>
// vector<T> matrix<T>::EigenValues(){
//     int MIN = -10;
//     int MAX = 10;
//     int sign = 0;
//     double STEP = 0.001;
//     double PRE = 0.01;

//     /*复制一份矩阵参与运算。*/
//     matrix<T> c_mat_alt = *this;
//     for(int real = MIN;real <= MAX;real += STEP)
//         for(int imag = MIN;imag <= MAX;imag += STEP)
//         {
//             /*以sign为0或1判断正的虚部之前是否需要输出加号。*/
//             sign = 0;
//             for(int i = 0;i < this->col;i ++){
//                 complex<double> temp (real,imag);
//                 c_mat_alt.data[i*this->col+i].operator-=(temp);
//             }
//                 // c_mat_alt.data[i*this->col+i] -= (real + imag * I);
//             /*如果行列式计算结果显示实部和虚部均小于给定精度（被认作为0），判定real + imag i是特征值。*/
//             if(fabs(creal(c_mat_alt.det())) <= PRE && fabs(cimag(c_mat_alt.det())) <= PRE)
//                 {
//                     /*调整输出格式，去掉不必要输出的部分。*/
//                     if(fabs(real) > PRE)
//                     {
//                         printf("%.3f ",real);
//                         sign = 1;
//                     }
//                     if(fabs(imag) > PRE)
//                     {
//                         if(imag > PRE)
//                         {
//                             if(sign == 1)
//                                 printf("+ ");
//                             printf("%.3fi",imag);
//                         }
//                         if(imag < - PRE)
//                             printf("%.3fi",imag);
//                     }
//                     if(fabs(real) < PRE && fabs(imag) < PRE)
//                         printf("0.000");
//                     printf("\n");
//                 }
//         }
// }

template <typename T>
vector<matrix<T>> matrix<T>::EigenVector()
{
    vector<T> eigenValue = this->EigenValues();

    matrix<T> c_ALU(this->row, this->col, this->data);

    for (int j = 0; j < eigenValue.size(); j++)
    {
        T cur = eigenValue.at(j);
        for (int i = 0; i < this->row + this->col; i++)
        {
            c_ALU.data[i] = this->data[i];
        }
        for (int i = 0; i < this->col; i++)
        {
            c_ALU.data[i * this->col + i] -= cur;
        }
        c_ALU = c_ALU.getRowEchelon();
        int rank = c_ALU.getRank();
        int numberOfValue = this->row - rank;
        vector<matrix<T>> values;

        // for(int i=0;i<numberOfValue;i++){
        //     values[i] = matrix(this->col,1);
        // }

        // int prow = 0;
        // int pcol = 0;
        // while (prow < c_ALU.col && pcol < c_ALU.row)
        // {
        //     if (c_ALU.data[prow * c_ALU.col + pcol] != 0)
        //     {
        //         prow++;
        //         pcol++;
        //     }
        //     else
        //     {
        //         T arry[c_ALU.col];
        //         for (int l = 0; l < c_ALU.row; l++)
        //         {
        //             arry[l] = (-1) * c_ALU.data[l * c_ALU.col + pcol];
        //         }
        //         arry[prow] = 1;
        //         values.push_back(matrix(this->col, 1, arry));
        //         pcol++;
        //     }
        // }

        int v = 0; //有n列无主元
        for (int jj = 0; jj < c_ALU.col; jj++)
        {
            for (int ii = jj - v; ii < c_ALU.row; ii++)
            { //对于第[j]列 本来是第[j]行 但因为跳过了v 我们的主元要从第[j-v]行开始找主元
                if (c_ALU.data[ii * this->col + jj] == 0)
                { //当前行主元列为0
                    if (ii == row - 1)
                    { //整列都是0 没有主元
                    T arry[c_ALU.col]{0};
                    for (int l = 0; l < c_ALU.row; l++)
                    {
                        if(c_ALU.data[l * c_ALU.col + jj]!=0)
                            arry[l] = (-1) * c_ALU.data[l * c_ALU.col + jj];
                    }
                    arry[jj] = 1;
                    values.push_back(matrix(this->col, 1, arry));
                        v++;
                        break;
                    }
                    else
                    {
                        continue;
                    }
                }
                else break;
            }
        }
        return values;
    }
}

// /*在这里定义和修改尝试的各种参数。*/
// #define MIN -10  /*尝试区间的下限*/
// #define MAX 10  /*尝试区间的上限*/
// #define STEP 0.001  /*尝试过程的步长*/
// #define PRE 0.01  /*判断为零的精度*/
// /*定义矩阵类型，其中C_MATRIX类型储存的元素均为复数。*/
// // typedef complex double C_MATRIX[20][20];
// typedef double MATRIX[20][20];
// /*定义复数类型的函数，用于计算行列式的值。*/
// complex double determinant(C_MATRIX c_mat,int order);
// complex double cofactor(C_MATRIX c_mat,int order,int r,int c);

// int main(int argc,char *argv[])
// {
//     double real,imag;
//     MATRIX mat;
//     C_MATRIX c_mat,c_mat_alt;
//     int order,i,j,sign = 0;

//     printf("输入矩阵的阶数:");
//     scanf("%d",&order);
//     printf("输入矩阵:\n");
//     for(i = 0;i < order;i ++)
//         for(j = 0;j < order;j ++)
//         {
//             scanf("%lf",&mat[i][j]);
//             /*将输入到mat中的元素赋值给存储复数类型元素的矩阵c_mat。*/
//             c_mat[i][j] = (complex double)mat[i][j];
//         }
//     /*试根求复特征值。*/
//     printf("特征值为:\n");
//     for(real = MIN;real <= MAX;real += STEP)
//         for(imag = MIN;imag <= MAX;imag += STEP)
//         {
//             /*以sign为0或1判断正的虚部之前是否需要输出加号。*/
//             sign = 0;
//             /*复制一份矩阵参与运算。*/
//             for(i = 0;i < order;i ++)
//                 for(j = 0;j < order;j ++)
//                     c_mat_alt[i][j] = c_mat[i][j];
//             for(i = 0;i < order;i ++)
//                 c_mat_alt[i][i] -= (real + imag * I);
//             /*如果行列式计算结果显示实部和虚部均小于给定精度（被认作为0），判定real + imag i是特征值。*/
//             if(fabs(creal(determinant(c_mat_alt,order))) <= PRE && fabs(cimag(determinant(c_mat_alt,order))) <= PRE)
//                 {
//                     /*调整输出格式，去掉不必要输出的部分。*/
//                     if(fabs(real) > PRE)
//                     {
//                         printf("%.3f ",real);
//                         sign = 1;
//                     }
//                     if(fabs(imag) > PRE)
//                     {
//                         if(imag > PRE)
//                         {
//                             if(sign == 1)
//                                 printf("+ ");
//                             printf("%.3fi",imag);
//                         }
//                         if(imag < - PRE)
//                             printf("%.3fi",imag);
//                     }
//                     if(fabs(real) < PRE && fabs(imag) < PRE)
//                         printf("0.000");
//                     printf("\n");
//                 }
//         }

//     return 0;
// }

template <typename T>
const T matrix<T>::getRank() const
{
    matrix<T> Matrix_RE = this->getRowEchelon();
    //得到简化阶梯矩阵来求秩
    // Matrix_RE.show();

    int rows = Matrix_RE.row;
    int cols = Matrix_RE.col;
    int rank = rows;
    for (int i = 0; i < rows; i++)
    {
        bool allzero = true;
        for (int j = 0; j < cols; j++)
        {
            if (Matrix_RE.data[i * this->col + j] != 0)
            {
                allzero = false;
                break;
            }
        }
        if (allzero)
        {
            rank--;
        }
    }
    return rank;
}

template <typename T>
matrix<T> matrix<T>::getRowEchelon() const
{
    matrix<T> Matrix_RE = matrix<T>(*this);
    int v = 0; //有n列无主元
    int factor = 1;
    for (int j = 0; j < Matrix_RE.col; j++)
    {
        for (int i = j - v; i < Matrix_RE.row; i++)
        { //对于第[j]列 本来是第[j]行 但因为跳过了v 我们的主元要从第[j-v]行开始找主元
            if (Matrix_RE.data[i * this->col + j] == 0)
            { //当前行主元列为0
                if (i == row - 1)
                { //整列都是0 没有主元
                    v++;
                    break;
                }
                else
                {
                    continue;
                }
            }
            else
            { //找到非零元 先交换到第j-v行 然后把其他都化为0
                Matrix_RE.RowExchange(i, j - v);
                // T pivot = Matrix_RE.vectors[j - v][j];//找到首元位置
                T pivot = Matrix_RE.data[(j - v) * this->col + j];
                for (int i_cur = 0; i_cur < row; i_cur++)
                {
                    T AtPivotCol = Matrix_RE.data[(i_cur) * this->col + j];
                    // T AtPivotCol = Matrix_RE.vectors[i_cur][j];//看第i_cur行的主元列是不是0 如果不是 每一项都乘以主元然后和主元*AtPivotCol相减成0
                    if (AtPivotCol != 0 && i_cur != j - v)
                    {
                        for (int j_cur = 0; j_cur < col; j_cur++)
                        {
                            Matrix_RE.data[(i_cur) * this->col + j_cur] *= pivot;
                            Matrix_RE.data[(j - v) * this->col + j_cur] *= AtPivotCol;
                            Matrix_RE.data[(i_cur) * this->col + j_cur] -= Matrix_RE.data[(j - v) * this->col + j_cur];
                            // Matrix_RE.vectors[i_cur][j_cur] /= pivot;
                            Matrix_RE.data[(j - v) * this->col + j_cur] /= AtPivotCol;
                        }
                    }
                }
                for (int o = 0; o < Matrix_RE.col; o++)
                {
                    Matrix_RE.data[(j - v) * this->col + o] /= pivot;
                }
                break;
            }
        }
    }

    return Matrix_RE;
}

template <typename T>
void matrix<T>::RowExchange(int a, int b)
{
    T temp[this->col];
    for (int i = 0; i < this->col; i++)
    {
        temp[i] = this->data[a * this->col + i];
    }
    for (int i = 0; i < this->col; i++)
    {
        this->data[a * this->col + i] = this->data[b * this->col + i];
    }
    for (int i = 0; i < this->col; i++)
    {
        this->data[b * this->col + i] = temp[i];
    }

    //     // vector<T> tempStore = vector<T>(this->col);
    //     tempStore = this->vectors[a];
    //     this->vectors[a] = this->vectors[b];
    //     this->vectors[b] = tempStore;
}

//传入图片地址，输出二维vector
template <typename T>
matrix<T> matrix<T>::decode(char* path)   //path为图片路径
{
    Mat img = imread(path);                // 将图片传入Mat容器中
//       显示原图片
//       namedWindow("old", WINDOW_NORMAL);
//       imshow("old", img);
//      waitKey(0);
    int w = img.cols * img.channels();     //可能为3通道，宽度要乘图片的通道数
    int h = img.rows;

    matrix<T> ar(h,w);
    // vector<vector<u_char> > array(h, vector<u_char>(w));      //初始化二维vector
    for (int i = 0; i < h; i++)
    {
        u_char* inData = img.ptr<u_char>(i);            //ptr为指向图片的行指针，参数i为行数
        for (int j = 0; j < w; j++)
        {
            ar.data[i*ar.col+j] = inData[j];
        }
    }
    return ar;
}

template <typename T>
cv::Mat matrix<T>::code(char* path)
{
    size_t h = array.size();                  
    size_t w = array[0].size();
    //初始化图片的像素长宽
    Mat img(h, (size_t)(w/3), CV_8UC3);           //保存为RGB，图像列数像素要除以3；
    for (size_t i = 0; i < h; i++)
    {
        u_char* outData = img.ptr<u_char>(i);
        for (size_t j = 0; j < w; j++)
        {
             outData[j] = this->data[i*this->col+j];
        }
    }
    namedWindow("new", WINDOW_NORMAL);
    imshow("new", img);
    waitKey(0);
    imwrite("目标地址\\save.jpg", img);
}


int main()
{
    // complex<double> a (3,0);
    // complex<double> b (-1,0);
    // complex<double> c (-1,0);
    // complex<double> d (3,0);
    // complex<double> arr[] = {3,-1,-1,3};
    // matrix<complex<double>> gg (2,2,arr);
    // double det = a.det();
    // cout<<det<<endl;
    // matrix<double> reverse = a.inverse();
    // reverse.display();
    // gg.values();

    double arr1[]{4, 1, 1, 4, 1, 1, 0, 0, 0};
    matrix<double> text(3, 3, arr1);
    text.getRowEchelon().display();
    cout << text.getRank() << endl;
    vector<matrix<double>> res = text.EigenVector();
    for (int i = 0; i < res.size(); i++)
    {
        res.at(i).display();
        cout << endl;
    }
}


