#include <stdio.h>
#include <stdlib.h>
#include <random>
#include "mex.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //declare variables
    mxArray *Nin, *Pin, *Rout;
    //int repeat;
    const mwSize *dims;
    double *n, *p, *r;
    int rows, cols;

    //associate inputs
    Nin = mxDuplicateArray(prhs[0]);
    Pin = mxDuplicateArray(prhs[1]);
    //repeat = mxGetScalar(prhs[2]);
    
    //figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    rows = (int)dims[0]; 
    cols = (int)dims[1];

    //associate outputs
    Rout = plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);

    //associate pointers
    n = mxGetPr(Nin);
    p = mxGetPr(Pin);
    r = mxGetPr(Rout);
    
    default_random_engine rng(random_device{}());
    for (int i = 0; i < rows * cols; i ++) {
        binomial_distribution<int> bino(n[i], p[i]);
        r[i] = bino(rng);
    }
    
    
//    if (repeat == 0) { // each row is unique
//        for (int i = 0; i < dim_i; i ++) {
//            for (int j = 0; j < dim_j; j ++) {
//                binomial_distribution<int> distribution(n[j * dim_i + i], p[j * dim_i + i]);
//                r[j * dim_i + i] = distribution(generator);
//            }
//        }
//    } 
//    else { // rows are identical
//        for (int j = 0; j < dim_j; j ++) {
//            binomial_distribution<int> distribution(n[j * dim_i], p[j * dim_i]);
//            for (int i = 0; i < dim_i; i ++) {
//                r[j * dim_i + i] = distribution(generator);
//            }
//        }
//    }
    
    return;
}