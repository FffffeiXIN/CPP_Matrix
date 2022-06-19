#ifndef CPP_PROJECT_TEST_CONVERSION_H
#define CPP_PROJECT_TEST_CONVERSION_H

#include <iostream>
#include<opencv2/opencv.hpp>
#include "matrix.h"


using namespace std;
using namespace cv;

template<typename T>
cv::Mat matrixToImag(matrix<T> mat) {
    size_t h = mat.get_row();
    size_t w = mat.get_col();
    //初始化图片的像素长宽
    Mat img(h, (size_t) (w / 3), CV_8UC3);     //保存为RGB，图像列数像素要除以3；
    for (size_t i = 0; i < h; i++) {
        uchar *outData = img.ptr<uchar>(i);
        for (size_t j = 0; j < w; j++) {
            outData[j] = (uchar) mat.getShape()[i * mat.get_col() + j];
        }
    }
    namedWindow("new", WINDOW_NORMAL);
    imshow("new", img);
    waitKey(0);
    imwrite("F:\\cpp_project_test\\save.jpg", img);
    return img;
}

matrix<int> imag_to_matrix(cv::Mat img) {
    int w = img.cols * img.channels();     //3通道，宽度要乘图片的通道数
    int h = img.rows;

    matrix<int> ar(h, w);
    for (int i = 0; i < h; i++) {
        uchar *inData = img.ptr<uchar>(i);    //ptr为指向图片的行指针，参数i为行数
        for (int j = 0; j < w; j++) {
            ar.getShape()[i * ar.get_col() + j] = (int) inData[j];
        }
    }
    return ar;
}


#endif //CPP_PROJECT_TEST_CONVERSION_H
