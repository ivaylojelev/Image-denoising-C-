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


function Lx = L( m )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Lx = zeros(m, m);
for i = 1:m-1
    Lx(i, i) = -1;
    Lx(i, i+1) = 1;

end

function U = PDHGMp ( F, rho, theta, sigma )
% That is the PDHGMp denoising alghotithm
% F is the original image, rho and theta are coefficients
% such that their product doesn't exceed 1/9
% and U is the clear image

K = 10000; %K is the maximum number of iterations of the optimizing algorithm
%eps = 1*e-10; % epsilon is a treshold for the conversion criterion

% Initializing U, P1, P2, P1_bar, P2_bar to be zero matrices of the size of
% F
[m, n] = size(F);
U = zeros(m, n);
P1 = zeros(m, n);
P2x = zeros(m, n);
P2y = zeros(m, n);
P1_bar = zeros(m, n);
P2x_bar = zeros(m, n);
P2y_bar = zeros(m, n);
%Set Lx and Ly
Lx = L(m);
Ly = L(n);

% Algorithm for optimizing;
% At each iteration U, P1, P2, P1_bar and P2_bar are updated until the
% algorithm converges with the approximation of the given threshold epsilon
for k = 1:K
    %First step of the algorithm is updating U
    DP2 = Lx' * P2x_bar + P2y_bar * Ly;
    X = U - (P1_bar + DP2) * rho * theta;
    U = min(max(0, X), 1);
    
    %Second step is projecting (u+p1) - point onto an N-dimensional sphere
    %centered at the input image F with a radius of R 
    R = sqrt(m*n*sigma);
    V1 = Ball_proj(U + P1, F, R);
    
    %Third step requires adding proximity mapping to our cost function
    %ROF and minimizing it by differentiating at zero
    Gx = P2x + Lx * U;
    Gy = P2y + U * Ly';
    Gr = sqrt(Gx.*Gx + Gy.*Gy);
    V2x = Gx - ((1/theta) * Gx)./max(Gr, 1/theta);
    V2y = Gy - ((1/theta) * Gy)./max(Gr, 1/theta);
    
    %Initializing the dual variables for the next iteration
    P1_prev = P1;
    P2_prev_x = P2x;
    P2_prev_y = P2y;
    P1 = P1 + U - V1;
    P2x = P2x + Lx * U - V2x;
    P2y = P2y + U * Ly' - V2y;
    
    P1_bar = P1 + ( P1 - P1_prev);
    P2x_bar = P2x + ( P2x - P2_prev_x);
    P2y_bar = P2y + ( P2y - P2_prev_y);
%    Conv = P1-P1_prev;
    
end
end
*/

using namespace cv;
using namespace std;


Mat<double> PDHGMp(Mat<double> F, double rho, double theta, double sigma)
{
    int K = 10000; // K is the maximum number of iterations of the optimizing algorithm
    
    int m = F.rows;
    int n = F.cols;
    U = Mat::zeros(m, n);
    P1 = Mat::zeros(m, n);
    P2x = Mat::zeros(m, n);
    P2y = Mat::zeros(m, n);
    P1_bar = Mat::zeros(m, n);
    P2x_bar = Mat::zeros(m, n);
    P2y_bar = Mat::zeros(m, n);
    
    // Set Lx and Ly
    Lx = Mat::zeros(m, m);
    for(int i = 0; i < m-2; i++)
    {
        Lx.ptr<double>(i)[i] = -1
        Lx.ptr<double>(i)[i+1] = 1;
    }

    Ly = Mat::zeros(n, n);
    for(int i = 0; i < n-2; i++)
    {
        Lx.ptr<double>(i)[i] = -1
        Lx.ptr<double>(i)[i+1] = 1;
    }

    for(int i = 0; i < K; i++)
    {
        // First step of the algorithm is updating U
        Mat<double> DP2 = transpose(Lx) * P2x_bar + P2y_bar * Ly;
        Mat<double> X = U - (P1_bar + DP2) * rho * theta;
        Mat<double> squeeze(X);
        
        //Second step is projecting (u+p1) - point onto an N-dimensional sphere
        // centered at the input image F with a radius of R 
        R = sqrt(m*n*sigma);
        V1 = Ball_proj(U + P1, F, R);
        
        // Third step requires adding proximity mapping to our cost function
        // ROF and minimizing it by differentiating at zero
        Gx = P2x + Lx * U;
        Gy = P2y + U * transpose(Ly);
        Gr = sqrt(sqr(Gx) + sqr(Gy));
        V2x = Gx - ((1/theta) * Gx) / max(Gr, 1/theta);
        V2y = Gy - ((1/theta) * Gy) / max(Gr, 1/theta);
        
        // Initializing the dual variables for the next iteration
        P1_prev = P1;
        P2_prev_x = P2x;
        P2_prev_y = P2y;
        P1 = P1 + U - V1;
        P2x = P2x + Lx * U - V2x;
        P2y = P2y + U * transpose(Ly) - V2y;
        
        P1_bar = P1 + ( P1 - P1_prev);
        P2x_bar = P2x + ( P2x - P2_prev_x);
        P2y_bar = P2y + ( P2y - P2_prev_y);
    }

    return U;
}