#ifndef MATRIX_SPARSEMATRIX_H
#define MATRIX_SPARSEMATRIX_H
// #include "complex.h"
#include <complex>
#include <utility>
#include "matrix.h"

#include "exception.h"

using namespace std;

/*
 * Store SparseMatrix
 */

struct position {
    int r;
    int c;

    position(int row, int col) : r(row), c(col) {};
};


template<typename T>
class SparseMatrix {
private:
    vector<T> values;
    vector<position> positions;
    int row;
    int col;


public:
    SparseMatrix(vector<T> vals, vector<position> pots, int a = 1, int b = 1) :
            values(vals), positions(move(pots)), row(a), col(b) {
//        cout << "The constructor of SpareMatrix is called" << endl;
        if (positions.size() != values.size()) {
            throw SpareSizeNotMatchException();
        }
        if ((double) values.size() / (row * col) > 0.2) {
            throw NotSpareException();
        }
    };


    SparseMatrix(matrix<T> &matrix) {
        this->row = matrix.get_row();
        this->col = matrix.get_col();
        long long cnt = 0;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                if (matrix.getShape()[i * matrix.get_col() + j]) {
                    cnt++;
                    positions.emplace_back(i, j);
                    values.push_back(matrix.getShape()[i * matrix.get_col() + j]);
                }
            }
        }
        //当矩阵没有达到“稀疏”的要求时，转换失败
        if ((double) cnt / (row * col) > 0.2) {
            throw NotSpareException();
        }
    };


    static matrix<T> convertToMatrix(SparseMatrix &sm) {
        matrix<T> matrix(sm.row, sm.col);
        for (int i = 0; i < sm.values.size(); i++) {
            matrix.getShape()[sm.positions[i].r * sm.col + sm.positions[i].c] = sm.values[i];
        }
        *matrix.get_ref_count() = 1;
        return matrix;
    }

    void dispaly() {
        cout << "The values are : ";
        for (int i = 0; i < values.size(); ++i) {
            cout << values[i] << " ";
        }
        cout << endl;
        cout << "The positions are : ";
        for (int i = 0; i < positions.size(); ++i) {
            cout << "(" << positions[i].r << "," << positions[i].c << ")" << " ";
        }
        cout << endl;
    }
};

#endif