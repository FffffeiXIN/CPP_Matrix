#include <iostream>
#include "matrix.cpp"

using namespace std;

void show(const vector<double> &v);

int main() {
    try {

        double array1[6] = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
        double array2[6] = {6.6, 5.5, 4.4, 3.3, 2.2, 1.1};
        int array3[6] = {1, 2, 3, 4, 5, 6};
        int array4[6] = {6, 5, 4, 3, 2, 1};
        //------test *
        {
            matrix<double> a(2, 3, array1);
            matrix<double> b(3, 2, array2);
            matrix<double> answer = a * b; //调用copy constructor
            answer.display();
            cout << "-----------------------------------" << endl;
        }
        //------test 2
        {
            matrix<double> a(3, 2, array1);
            matrix<double> b(3, 2, array2);
            (a + b).display();
            (a - b).display();
            (a * 2).display();
            (2.0 * a).display();
            (a / 2.0).display();
            cout << a.dot(b) << endl;
            cout << "-----------------------------------" << endl;
        }
        {
            matrix<int> a(3, 2, array3);
            (a / 2).display();
            cout << "-----------------------------------" << endl;
        }

        //------test 3 assign copy and get max or min or average
        {
            matrix<double> b(3, 2, array1);
            matrix<double> a = b;
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
            v = a.average("row");
            show(v);
            v = a.average("col");
            show(v);
            v = a.average();
            show(v);
            cout << "-----------------------------------" << endl;
        }

        //------test 4
        {
            matrix<double> a(3, 2, array1);
            a.transpose().display();
            a.reshape(6, 1).display();

            a.slicing(0, 2, 0, 2).display();
            a.slicing(1, 1).display();
            cout << "-----------------------------------" << endl;
        }
        {
            double arr[] = {1.0, 2, 3, 2, 2, 1, 3, 4, 3};
            matrix<double> a(3, 3, arr);

            double det = a.det();
            cout << det << endl;

            matrix<double> reverse = a.inverse();
            reverse.display();
            cout << "-----------------------------------" << endl;
        }
        {
            int array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
            matrix<int> a(3, 3, array);
            cout << a.trace() << endl;
            cout << "-----------------------------------" << endl;
        }
        {
            matrix<double> a(3, 2, array1);
            vector<double> v = {3.6, 7.3};
            a.vector_multiplication(v).display();
            cout << "-----------------------------------" << endl;
        }
        {
            matrix<double> a(3, 2, array1);
            vector<double> v = {3.6, 7.3, 2.2};
            matrix<double>::vector_multiplication(v, a).display();
            cout << "-----------------------------------" << endl;
        }
    }
    catch (exception &e) {
        cout << e.what() << endl;
    }
    return 0;
}

void show(const vector<double> &v) {
    for (double i: v) {
        cout << i << " ";
    }
    cout << endl;
}
