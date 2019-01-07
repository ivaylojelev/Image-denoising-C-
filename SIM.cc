#include <vector>
#include <stdio.h>
#include <string.h>
#include "mean_var.h"
#include "Get_Intersection.cpp"
#include "Get_Sphere.cpp"
#include "Ball_proj.cpp"
#include "Ball_proj_SIM.cpp"
#include "PDHGMp.cpp"
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

Dim = size(F(:),1);
R1 = sqrt((sigma)*Dim);
[F1, R] = Get_Sphere(F, sigma);
sigma1 = var(F1(:)-F(:));
sigma_curr = sigma+sigma1;
for i=2:Sph_n
    [FSC, RSC, sigma_curr] = Get_Sphere(FSC, sigma_curr);
    [InterS, InterR, B] = Get_Intersection(FSC, RSC, B, InterS, InterR);
   end
S = reshape(InterS(Sph_n,:,:),m,n);

U = Ball_proj_SIM(U_P, FSC, RSC, B);
end

In this version, instead of taking the noisy image as an argument and the clear image as a return value,
we take an url string for both images. Then the algorithm reads the image with OpenCV functionality
and saves the clear one

*/

using namespace cv;
using namespace std;

int main()
{
    string urlF, urlU;
    double sigma;
    
    cout<<"Input the url of the noisy image:"<<endl;
    cin.getline(urlF);
    cout<<"Input the url for saving the clear image:"<<endl;
    cin.getline(urlU);
    
    cout<<"Input the noise variance sigma:"<<endl;
    cin>>sigma;

    Mat<double> F = imread(urlF, CV_LOAD_IMAGE_GRAYSCALE);
    F /= 255;

    Mat<double> F1 = zeros(F.size()), FSC = zeros(F.size()), InterS = zeros(F.size());
    vector<Mat<double> B = new vector<Mat<double>>();
    double dim = F.rows * F.cols;
    double sigma1, R, InterR;

    F1 = Get_Sphere(&R, &sigma1, F, sigma);
    sigma1 = Var(F1 - F);

    sigma_curr = sigma + sigma1;
    R = sqrt(dim*sigma1);
    InterR = R;

    for(int i = 2; i <= Sph_n; i++) {
        FSC = Get_Sphere(&RSC, &sigma_curr, FSC, sigma_curr);
        InterS = Get_Intersection(&InterR, &B, B, InterS, InterR);
    }
    Mat<double> U_P = PDHGMp(F, 0.2, 0.5, sigma);
    Mat<double> U = Ball_proj_SIM(U_P, FSC, RSC, B);

    return imwrite(urlU, U.data);
}