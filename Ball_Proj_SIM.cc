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


[k,m,n] = size(B);
proj = P;
for i = 1:k 
    vec = proj-S;
    proj = proj - sum(dot(vec,reshape(B(i,:,:),m,n))) * reshape(B(i,:,:),m,n);
end
vec = proj-S;
M = Ball_SIM(proj, S, R)
*/

using namespace cv;
using namespace std;


Mat<double> Ball_proj(Mat<double> P, Mat<double> P, double R, vector<Mat<double>> B)
{
    int k = B.size();
    Mat<double> proj = P;
    Mat<double> vec = zeros(P.size());
    for(int i = 0; i < k; i++)
    {
        vec = proj - S;
        proj = proj - sum(dot(vec, B[i]) * B[i])
    }
    vec = proj - S;
    return Ball_SIM(proj, S, R);
}