#ifndef MEANVAR
#define MEANVAR
#endif

#include <vector>
#include <stdio.h>
#include <math.h>
#include <opencv2/opencv.hpp>
double Mean(Mat<double>m)
{
    double sum = 0;
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.cols; j++) sum += m.ptr<double>(i)[j];
    }
    return sum / (m.rows * m.cols);
}
double Var(Mat<double>m)
{
    double sum = 0;
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.cols; j++) sum += m.ptr<double>(i)[j] * m.ptr<double>(i)[j];
    }
    return sqrt(sum / (m.rows * m.cols));
}

double sum(Mat<double>m)
{
    double sum = 0;
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.cols; j++) sum += m.ptr<double>(i)[j] * m.ptr<double>(i)[j];
    }
    return sum;
}

double dot(Mat<double> a, Mat<double> b)
{
    return sum(a*b);
}


Mat<double> sqr(Mat<double> m)
{
    Mat<double> ret = m;
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.cols; j++) ret.ptr<double>(i)[j] = m.ptr<double>(i)[j] * m.ptr<double>(i)[j];
    }
    return ret;
}

Mat<double> sqeeze(Mat<double> m)
{
    for (int i = 0; i < m.rows; i++)
    {
        for (int j = 0; j < m.cols; j++) 
        {
            if(m.ptr<double>(i)[j] < 0) m.ptr<double>(i)[j] = 0;
            if(m.ptr<double>(i)[j] > 1) m.ptr<double>(i)[j] = 1;
        }
    }
    return m;
}
