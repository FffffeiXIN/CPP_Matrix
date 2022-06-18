#ifndef CPP_PROJECT_TEST_EXCEPTION_H
#define CPP_PROJECT_TEST_EXCEPTION_H

#include <iostream>

using namespace std;

//两个矩阵或一个矩阵一个向量做运算时大小不匹配
class SizeMismatchException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The size of two matrices (vector) are not matching.\n";
    }
};

//除法不能除以0
class ZeroDivideException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The divisor is 0.\n";
    }
};

//Reshape的大小不合法
class InvalidReshapeSizeException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The new shape is illegal.\n";
    }
};

//该矩阵不是一个方阵
class NonSquareException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The matrix isn't a square matrix.\n";
    }
};

//矩阵不可逆
class IrreversibleException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The matrix is not invertible.\n";
    }
};

//矩阵不是一个3维的向量
class Not3DVectorException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The matrix is not a 3-dimensional vector .\n";
    }
};

//输入一些类型进而做操作时，类型不合法
class TypeInvalidException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The input type is invalid and cannot be operated on.\n";
    }
};

//稀疏矩阵的position和values的大小不相等，无法构建稀疏矩阵
class SpareSizeNotMatchException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The size of values and positions not match, the constructor can't work.\n";
    }
};

//该矩阵不是一个稀疏矩阵
class NotSpareException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The matrix is not a spare matrix.\n";
    }
};

#endif // CPP_PROJECT_TEST_EXCEPTION_H
