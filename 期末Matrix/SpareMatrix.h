#ifndef MATRIX_SPARSEMATRIX_H
#define MATRIX_SPARSEMATRIX_H
// #include "complex.h"
#include <complex>
#include "matrix.h"

/*
 * Store SparseMatrix
 */

struct position
{
    int r;
    int c;
    position(int row,int col):r(row),c(col){};
};


template <typename T>
class SparseMatrix{
private:
    vector<T> values;
    vector<position> positions;
    int row;
    int col;
    
    
public:
    SparseMatrix(vector<T> vals,vector<position> pots,int a=1,int b=1):values(vals),positions(pots),row(a),col(b){};
    SparseMatrix(matrix<T> & matrix){
        this->row=matrix.row;
        this->col=matrix.col;
        long long cnt=0;
        for(int i=0;i<row;i++){
            for (int j=0;j<col;j++){
                if(matrix.data[i*matrix.col+j]){
                    cnt++;
                    positions.emplace_back(i,j);
                    values.push_back(matrix.data[i*matrix.col+j]);
                }
            }
        }
        if((double)cnt/(row*col)>0.2){
            cerr<<"The matrix is not sparse enough.\n"
                  "The size is nearly the same.\n"
                  " So the suggestion is that use normal matrix to store it."<<endl;
        }
    };

    static matrix<T> convertToMatrix(SparseMatrix & sparseMatrix){
        matrix<T> matrix(sparseMatrix.row,sparseMatrix.col);
        for(int i=0;i<sparseMatrix.values.size();i++){
            matrix.vectors[sparseMatrix.positions[i].x*sparseMatrix.col+sparseMatrix.positions[i].y]
                    =sparseMatrix.values[i];
        }
        return matrix;
    }
};
#endif 
