# Image-denoising-C-
For the last year and a half I have been developing an innovative algorithm for extracting digital noise from images. The method has better results than present denoising algorithms and is based on Sphere-constrained total variation. We needed parallel processing for this so we used MatLAB scripts but we wanted to further improve the speed. Therefore, we translated the relevant pieces of code in C++. Translating a script from an interpreted, non-object-oriented language to compiled object-oriented language such as C++ requires fluency in both languages and poses a lot of problems mostly with memory allocation which I managed to handle well. We also needed to preserve the dynamic processing so we had to use the matrices operation functionality from OpenCV and implement it. 

The nature of the algorithm is very technical but if you want to learn more you can read the research paper (publication pending) I wrote on it: https://drive.google.com/file/d/1Hhc3AaXgMZsa7ZKHJScoGGPqeKNV-Ifa/view?usp=sharing



P.S. For developing this algorithm, I got awarded with the Bulgarian Presidential award "John Atanasov" (Oct 2018) and recieved two special awards in the European Union Contest for Young Scientists (Dublin, Sep 2018).
