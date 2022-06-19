#include <iostream>
#include "matrix.cpp"
#include "SpareMatrix.h"
#include <vector>

using namespace std;

template<typename T>
void show(const vector<T> &v);

int main() {
    try {
//        {
//            //-----test for SpareMatrix-----
//            vector<double> value = {1, 5, 10, 14};
//            position a(4, 3);
//            position b(9, 9);
//            position c(1, 3);
//            position d(5, 8);
//            vector<position> positions = {a, b, c, d};
//            SparseMatrix<double> sm(value, positions, 10, 10);
//            matrix<double> m = SparseMatrix<double>::convertToMatrix(sm);
//            m.display();
//            SparseMatrix<double> sm1(m);
//            cout << endl;
//            sm1.dispaly();
//        }
//        {
//            // test for value and position of SpareMatrix exception
//            vector<double> value = {1, 5, 10, 14};
//            position a(4, 3);
//            position b(9, 9);
//            position c(1, 3);
//            vector<position> positions1 = {a, b, c};
//            SparseMatrix<double> sm2(value, positions1, 10, 10);
//        }
//        cout << "----------------" << endl;
//        {
//            //-----test for constructor and destructor-----
//            double array1[6] = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
//            matrix<double> a(2, 3, array1);
//            a.display();
//
//            matrix<double> b(a);
//            b.display();
//            matrix<double> c(2, 3);
//            c = a;
//            c.display();
//            cout << a.getdata() << endl;
//            cout << b.getdata() << endl;
//            cout << c.getdata() << endl;
//            cout << a.get_ref_count() << endl;
//            cout << b.get_ref_count() << endl;
//            cout << c.get_ref_count() << endl;
//        }
//        cout << "----------------" << endl;
//        {
//            //-----test for NotSpareMatrix exception-----
//            vector<double> value = {1, 5, 10, 14};
//            position a(4, 3);
//            position b(9, 9);
//            position c(1, 3);
//            position d(5, 8);
//            vector<position> positions = {a, b, c, d};
//            SparseMatrix<double> sm(value, positions, 3, 3);
//        }
//        cout << "----------------" << endl;
//        {
//            //-----test for NotSpareMatrix during converting exception-----
//            int array[9] = {0, 0, 0, 1, 3, 0, 0, 4, 5};
//            matrix<int> m(3, 3, array);
//            SparseMatrix<int> sm(m);
//            sm.dispaly();
//            m.display();
//        }
        {
            // test for matrix-matrix multiplication
            double array1[6] = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
            double array2[6] = {6.6, 5.5, 4.4, 3.3, 2.2, 1.1};
            matrix<double> a(2, 3, array1);
            matrix<double> b(3, 2, array2);
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "---------- b is ----------" << endl;
            b.display();
            matrix<double> answer = a * b; //调用copy constructor
            cout << "---------- a * b is ----------" << endl;
            answer.display();
        }
        {
            double array1[6] = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
            double array2[6] = {6.6, 5.5, 4.4, 3.3, 2.2, 1.1};
            matrix<double> a(3, 2, array1);
            matrix<double> b(3, 2, array2);
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "---------- b is ----------" << endl;
            b.display();
            cout << "---------- a + b is ----------" << endl;
            (a + b).display();
            cout << "---------- a - b is ----------" << endl;
            (a - b).display();
            cout << "---------- a * 2 is ----------" << endl;
            (a * 2).display();
            cout << "---------- 2.0 * a is ----------" << endl;
            (2.0 * a).display();
            cout << "---------- a / 2.0 is ----------" << endl;
            (a / 2.0).display();
            cout << "---------- a.element_wise_multiplication(b) ----------" << endl;
            a.element_wise_multiplication(b).display();
            cout << "---------- the transpose of a is ----------" << endl;
            a.transpose().display();
            cout << "-----------------------------------" << endl;
        }
        {
            int array1[3] = {1, 2, 3};
            int array2[6] = {6, 5, 4, 3, 2, 1};
            matrix<int> a(3, 1, array1);
            matrix<int> b(3, 2, array2);
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "---------- b is ----------" << endl;
            b.display();
            cout << "---------- the dot of a and b is ----------" << endl;
            a.dot(b).display();
        }
        {
            int array1[3] = {0, 0, 1};
            int array2[3] = {0, 1, 0};
            matrix<int> a(1, 3, array1);
            matrix<int> b(1, 3, array2);
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "---------- b is ----------" << endl;
            b.display();
            cout << "---------- the cross multiplication of a and b is ----------" << endl;
            a.cross(b).display();
        }
        {
            int array[6] = {6, 5, 4, 3, 2, 1};
            matrix<int> a(3, 2, array);
            vector<int> v1 = {4, 4, 4};
            vector<int> v2 = {3, 3};
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "---------- matrix-vector multiplication of v1 * a is ----------" << endl;
            vector_multiplication(v1, a).display();
            cout << "---------- matrix-vector multiplication of a * v2 is ----------" << endl;
            a.vector_multiplication(v2).display();
        }
        {
            //------test for get max or min or average
            double array1[6] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
            matrix<double> a(3, 2, array1);
            vector<double> v;
            v = a.maxValue();
            show(v);
            v = a.minValue();
            show(v);
            v = a.maxValue("row");
            show(v);
            v = a.maxValue("col");
            show(v);
            v = a.minValue("row");
            show(v);
            v = a.minValue("col");
            show(v);
            v = a.sum();
            show(v);
            v = a.average("row");
            show(v);
            v = a.average("col");
            show(v);
            v = a.average();
            show(v);
        }
        {
            int array[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
            matrix<int> a(3, 3, array);
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "The trace of a is : " << a.trace() << endl;
        }
        {
            double array1[6] = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
            matrix<double> a(3, 2, array1);
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "---------- reshape a ----------" << endl;
            a.reshape(6, 1).display();
            cout << "---- Take the first 2 rows, the first 2 columns of a ----" << endl;
            a.slicing(0, 2, 0, 2).display();
            cout << "---- Take from 1st row, 1st col of a ----" << endl;
            a.slicing(1, 1).display();
        }
        {
            // Test for convolution
            int arr[] = {17, 24, 1, 8, 15,
                         23, 5, 7, 14, 16,
                         4, 6, 13, 20, 22,
                         10, 12, 19, 21, 3,
                         11, 18, 25, 2, 9};
            int kernel[] = {1, 3, 1, 0, 5, 0, 2, 1, 2};
            matrix<int> a(5, 5, arr);
            matrix<int> b(3, 3, kernel);
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "---------- b is ----------" << endl;
            b.display();
            cout << "---------- valid convolution----------" << endl;
            a.convolution(b, "valid").display();
            cout << "---------- same convolution----------" << endl;
            a.convolution(b, "same").display();
            cout << "---------- full convolution----------" << endl;
            a.convolution(b, "full").display();
        }
        {
            double arr[] = {1.0, 2.0, 3.0, 2.0, 2.0, 1.0, 3.0, 4.0, 3.0};
            matrix<double> a(3, 3, arr);
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "---- The determinant of a is: ";
            double det = a.det();
            cout << det << endl;
            cout << "---- The reverse of a is ----";
            matrix<double> reverse = a.inverse();
            reverse.display();
        }
        {
            // cv::Mat img = imread("F:\\cpp_project_test\\sky.png");
            // matrix<int> toMatrix = imag_to_matrix(img);
            // cout << *toMatrix.get_ref_count();
            // toMatrix.display();
            // matrixToImag(toMatrix);
        }
        {
            double arr1[]{3, 2, 4, 2, 0, 2, 4, 2, 3};
            matrix<double> a(3, 3, arr1);
            cout << "---------- a is ----------" << endl;
            a.display();
            cout << "---- The determinant of a is: ";
            double det = a.det();
            cout << det << endl;
            matrix<double> reverse = a.inverse();
            cout << "--- The eigenValues of a is ---" << endl;
            show(a.EigenValues());
            vector<matrix<double>> res = a.EigenVector();
            cout << "---------- The eigenVector of a are ----------" << endl;
            for (auto &re: res) {
                re.display();
                cout << endl;
            }
        }
        {
            complex<int> a(3);
            complex<int> b(2, -1);
            complex<int> c(2, 1);
            complex<int> d(1);
            complex<int> e(5, 8);
            complex<int> f(9, 0);
            complex<int> data_temp[6] = {a, b, c, d, e, f};
            matrix<complex<int>> m(3, 2, data_temp);
            cout << "---------- m is ----------" << endl;
            m.display_complexMatrix();
            cout << "---------- the conjugation of m is ----------" << endl;
            m.conjugation().display_complexMatrix();
        }
        {
            // test for matrix-matrix multiplication
            complex<double> a(3.5);
            complex<double> b(2, -1);
            complex<double> c(2.5, 1);
            complex<double> d(1);
            complex<double> e(5, 8.5);
            complex<double> f(9.5, 0);
            complex<double> data_temp1[6] = {a, b, c, d, e, f};
            complex<double> data_temp2[6] = {f, e, d, c, a, b};
            matrix<complex<double>> m(3, 2, data_temp1);
            matrix<complex<double>> n(2, 3, data_temp2);
            cout << "---------- m is ----------" << endl;
            m.display_complexMatrix();
            cout << "---------- n is ----------" << endl;
            n.display_complexMatrix();
            cout << "---------- m * n is ----------" << endl;
            (m * n).display_complexMatrix();
        }
        {
            complex<double> a(3.5);
            complex<double> b(2, -1);
            complex<double> c(2.5, 1);
            complex<double> d(1);
            complex<double> e(5, 8.5);
            complex<double> f(9.5, 0);
            complex<double> data_temp1[6] = {a, b, c, d};
            complex<double> data_temp2[6] = {f, e, d, c};
            matrix<complex<double>> m(2, 2, data_temp1);
            matrix<complex<double>> n(2, 2, data_temp2);
            cout << "---------- m is ----------" << endl;
            m.display_complexMatrix();
            cout << "---------- n is ----------" << endl;
            n.display_complexMatrix();
            cout << "---------- m + n is ----------" << endl;
            (m + n).display_complexMatrix();
            cout << "---------- m - n is ----------" << endl;
            (m - n).display_complexMatrix();
            cout << "---------- m * 2 is ----------" << endl;
            (m * 2).display_complexMatrix();
            cout << "---------- 2.0 * m is ----------" << endl;
            (2 * m).display_complexMatrix();
            cout << "---------- m / 2.0 is ----------" << endl;
            (m / 2.0).display_complexMatrix();
            cout << "---------- m.element_wise_multiplication(m) ----------" << endl;
            m.element_wise_multiplication(n).display_complexMatrix();
            cout << "---------- the transpose of m is ----------" << endl;
            m.transpose().display_complexMatrix();
        }
        {

            complex<double> a(3.5);
            complex<double> b(2, -1);
            complex<double> c(2.5, 1);
            complex<double> d(1);
            complex<double> e(5, 8.5);
            complex<double> f(9.5, 0);
            complex<double> data_temp1[6] = {a, b, c};
            complex<double> data_temp2[6] = {f, e, d, c, a, b};
            matrix<complex<double>> m(3, 1, data_temp1);
            matrix<complex<double>> n(3, 2, data_temp2);
            cout << "---------- m is ----------" << endl;
            m.display_complexMatrix();
            cout << "---------- n is ----------" << endl;
            n.display_complexMatrix();
            cout << "---------- the dot of m and n is ----------" << endl;
            (m.dot(n)).display_complexMatrix();
        }
        {
            // test for matrix-matrix multiplication
            complex<double> a(1, 1);
            complex<double> b(0, 0);
            complex<double> data_temp1[6] = {a, b, b};
            complex<double> data_temp2[6] = {b, a, b};
            matrix<complex<double>> m(1, 3, data_temp1);
            matrix<complex<double>> n(1, 3, data_temp2);
            cout << "---------- m is ----------" << endl;
            m.display_complexMatrix();
            cout << "---------- n is ----------" << endl;
            n.display_complexMatrix();
            cout << "---------- the cross of m and n is ----------" << endl;
            (m.cross(n)).display_complexMatrix();
        }
        {
            //------test for get max or min or average
            complex<double> a(3.5);
            complex<double> b(2, -1);
            complex<double> c(2.5, 1);
            complex<double> d(1);
            complex<double> e(5, 8.5);
            complex<double> f(9.5, 0);
            complex<double> data_temp1[6] = {a, b, c, d, e, f};
            matrix<complex<double>> m(3, 2, data_temp1);
            cout << "----- m is ----" << endl;
            m.display_complexMatrix();
            cout << "---------------" << endl;
            vector<complex<double>> v;
            v = m.sum();
            show(v);
            v = m.average("row");
            show(v);
            v = m.average("col");
            show(v);
            v = m.average();
            show(v);
        }
        {
            complex<double> a(3.5);
            complex<double> b(2, -1);
            complex<double> c(2.5, 1);
            complex<double> d(1);
            complex<double> data_temp1[4] = {a, b, c, d};
            matrix<complex<double>> m(2, 2, data_temp1);
            cout << "---------- m is ----------" << endl;
            m.display_complexMatrix();
            cout << "--------------------------" << endl;
            cout << "The trace of m is : ";
            display_complex(m.trace());
            cout << endl;
        }
        {
            complex<double> a(3.5);
            complex<double> b(2, -1);
            complex<double> c(2.5, 1);
            complex<double> d(1);
            complex<double> e(5, 8.5);
            complex<double> f(9.5, 0);
            complex<double> data_temp1[6] = {a, b, c, d, e, f};
            matrix<complex<double>> m(3, 2, data_temp1);
            cout << "---------- m is ----------" << endl;
            m.display_complexMatrix();
            cout << "---------- reshape m ----------" << endl;
            m.reshape(6, 1).display_complexMatrix();
            cout << "---- Take the first 2 rows, the first 2 columns of m ----" << endl;
            m.slicing(0, 2, 0, 2).display_complexMatrix();
            cout << "---- Take from 1st row, 1st col of m ----" << endl;
            m.slicing(1, 1).display_complexMatrix();
        }
        {
            complex<double> a(3, 6);
            complex<double> b(-1, 2);
            complex<double> c(-1, 1);
            complex<double> d(3, 0);
            complex<double> arr[] = {a, b, c, d};
            double arri[] = {3, -1, -1, 3};
            matrix<complex<double>> gg(2, 2, arr);
            matrix<double> hh(2, 2, arri);
            complex<double> det1 = gg.det_complex();
            double det2 = hh.det();
            cout << det1 << endl;
            cout << det2 << endl;
            gg.inverse_complex();
            hh.inverse();
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