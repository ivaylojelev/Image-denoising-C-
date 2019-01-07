#include <vector>
#include <stdio.h>
#include "mean_var.h"
#include <opencv2/opencv.hpp>

/* 
This is the original MatLAB script which uses adds noise to the noisy image to construct a second sphere
This code is emulated in C++ for speed purposes

We want to keep the parallel processing features, so we use the cv::Mat object for operating with matrices
and some other functionality from OpenCV

This serves to demonstrate how a code can be emulated to C++ from a very different language (not even
object-oriented) and how different techniques and functionality have to be applied

The algorithm is complicated so if you want to delve into its technicalities I suggest you read the full
research paper this was based on:https://drive.google.com/file/d/1Hhc3AaXgMZsa7ZKHJScoGGPqeKNV-Ifa/view?usp=sharing


function [ S, R, B ] = Get_Intersection( C, R1, B_prev, D, R2 )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
len = norm(C-D, 'fro');
B = B_prev;
B(end+1,:,:) = (D - C)/len;
a = (R2^2 - R1^2 + len^2)/(2*len^2);
S = D - a*(D-C);
len1 = norm(D-S, 'fro');
R = sqrt(R2^2-len1^2);
end

*/

using namespace cv;
using namespace std;

Mat<double> Get_Intersection(double &R, vector<Mat<double>> &B, vector<Mat<double>> B_prev,  Mat<double> C, Mat<double> D, double R1, double R2)
{
    double len = Var(C-D);
    vector<Mat<double>> B = B_prev;
    Mat<double> temp = new Mat<double>(C.size());
    temp = (C-D);

    for (int i = 0; i < temp.rows; i++)
    {
        for (int j = 0; j < temp.cols; j++) temp.ptr<double>(i)[j] /= len;
    }
    B->push_back(temp);

    double a = (R2*R2 - R1*R1 + len*len)/(2*len*len);
    Mat<double> S = new Mat<double>(C.size());

    for (int i = 0; i < temp.rows; i++)
    {
        for (int j = 0; j < temp.cols; j++) temp.ptr<double>(i)[j] *= (len*a);
    }
    S = D - temp;
    
    len1 = Var(D-S);
    R = sqrt(R2*R2 - len1*len1);

    return S;
}