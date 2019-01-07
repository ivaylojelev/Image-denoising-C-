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


function M = Ball_proj( P, S, R )

vec = P-S;
len = norm(vec, 'fro');
if (len < R) 
    M = P;
else M = (R/len)*vec + S;

end

*/

using namespace cv;
using namespace std;


Mat<double> Ball_proj(Mat<double> P, Mat<double> P, double R)
{
    Mat<double> M = zeros(P.size());
    Mat<double> vec = P-S;
    len = Var(vec);

    if(len < R) M = P;
    else M = (R/len)*vec + S;

    return M;
}