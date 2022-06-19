#include <iostream>
#include "matrix.cpp"
#include "SpareMatrix.h"
#include "conversion.h"

using namespace std;

template<typename T>
void show(const vector<T> &v);

int main() {


    try {
        {
            complex<int> a(3);
            complex<int> b(2, -1);
            complex<int> c(2, 1);
            complex<int> d(1);
            complex<int> e(5, 8);
            complex<int> f(9, 0);
            complex<int> data_temp1[6] = {a, b, c, d};
            complex<int> data_temp2[6] = {f, e, d, c};
            matrix<complex<int>> m(2, 2, data_temp1);
            matrix<complex<int>> n(2, 2, data_temp2);
            cout << "---------- m is ----------" << endl;
            m.display_complexMatrix();
            cout << "---------- n is ----------" << endl;
            n.display_complexMatrix();
            cout << "---------- a + b is ----------" << endl;
            (m + n).display_complexMatrix();
            cout << "---------- a - b is ----------" << endl;
            (m - n).display_complexMatrix();
//            cout << "---------- a * 2 is ----------" << endl;
//            (m * 2).display_complexMatrix();
//            cout << "---------- 2.0 * a is ----------" << endl;
//            (2 * m).display_complexMatrix();
//            cout << "---------- a / 2.0 is ----------" << endl;
//            (m / 2.0).display_complexMatrix();
            cout << "---------- a.element_wise_multiplication(b) ----------" << endl;
            m.element_wise_multiplication(n).display_complexMatrix();
            cout << "---------- the transpose of a is ----------" << endl;
            m.transpose().display_complexMatrix();

        }
    }
    catch (exception &e) {
        cout << e.what() << endl;
    }
    return 0;
}

template<typename T>
void show(const vector<T> &v) {
    for (T i: v) {
        cout << i << " ";
    }
    cout << endl;
}
