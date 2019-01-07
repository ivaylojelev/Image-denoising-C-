#include <vector>
#include <stdio.h>
#include <math.h>
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


function [ F_next, R, sigma_bar ] = Get_Sphere( F_n, sigma_prev )
dim = size(F_n(:),1);
sigma1 = 0.01;
U = zeros(size(F_n))+0.5;
mu = 0.01;
N1 = imnoise(U, 'gaussian', 0, sigma1) - U + mu;
N1 = N1 - mean(N1(:));
F_next = F_n + N1;
sigma = var(N1(:));
sigma_bar = sigma_prev + sigma;
%(-b-sqrt(D))/2*a;
%Norm = norm(F1 - U, 'fro');
%Norm
R = sqrt((sigma_bar) * dim);

end
*/

using namespace cv;
using namespace std;

Mat Get_Sphere(double &R, double &sigma_bar, Mat<double> F_n, double sigma_prev)
{
    int dim = F_n.rows * F_n.cols;
    double sigma1 = 0.01;
    U = Mat::zeros(cv::Size(F_n.rows, F_n.cols);
    for (int i = 0; i < U.rows; i++)
    {
        for (int j = 0; j < U.cols; j++) U.ptr<double>(i)[j] = 0.5;
    }
    double mu = 0.01;
    Mat<double> noise(F_n.size(),double);
    Mat N1 = randn(noise, mu, sigma1) - U;
    N1 -= Mean(N1);
    
    Mat<double> F_next = F_n + N1;
    double sigma = Var(N1);
    
    sigma_bar = sigma + sigma_prev;
    double Norm = Var(F1 - U);
    cout<<Norm<<endl;

    R = sqrt((sigma_bar) * dim);
    return F_next;
}
