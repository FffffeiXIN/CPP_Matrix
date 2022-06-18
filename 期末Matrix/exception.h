#ifndef CPP_PROJECT_TEST_EXCEPTION_H
#define CPP_PROJECT_TEST_EXCEPTION_H

#include <iostream>

using namespace std;

class SizeMismatchException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The size of two matrices (vector) are not matching.\n";
    }
};

class ZeroDivideException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The divisor is 0.\n";
    }
};

class InvalidReshapeSizeException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The new shape is illegal.\n";
    }
};

class NonSquareException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The matrix isn't a square matrix.\n";
    }
};

class IrreversibleException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The matrix is not invertible.\n";
    }
};

class Not3DVectorException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The matrix is not a 3-dimensional vector .\n";
    }
};

class TypeInvalidException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The input type is invalid and cannot be operated on.\n";
    }
};

class SpareSizeNotMatchException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The size of values and positions not match, the constructor can't work.\n";
    }
};

class NotSpareException : public exception {
public:
    const char *what() const noexcept override {
        return "Exception: The matrix is not a spare matrix.\n";
    }
};

#endif // CPP_PROJECT_TEST_EXCEPTION_H
