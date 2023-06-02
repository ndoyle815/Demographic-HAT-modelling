#include <stdio.h>
#include <stdlib.h>
#include <random>
#include "mex.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //declare variables
    mxArray *Nin, *Pin, *Rhoin, *Rout;
    const mwSize *dims;
    double *n, *p, *rho, *r;
    int rows, cols;

    //associate inputs
    Nin = mxDuplicateArray(prhs[0]);
    Pin = mxDuplicateArray(prhs[1]);
    Rhoin = mxDuplicateArray(prhs[2]);
    
    //figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    rows = (int)dims[0]; 
    cols = (int)dims[1];

    //associate outputs
    Rout = plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);

    //associate pointers
    n = mxGetPr(Nin);
    p = mxGetPr(Pin);
    rho = mxGetPr(Rhoin);
    r = mxGetPr(Rout);
    
    
    default_random_engine rng(random_device{}());
    for (int i = 0; i < rows * cols; i ++) {
        if (n[i] == 0 || p[i] == 0) {
            r[i] = 0;
        } 
        else if (rho[i] == 0) {
            binomial_distribution<int> bino(n[i], p[i]);
            r[i] = bino(rng);
        }
        else {
            gamma_distribution<double> gamma1(p[i] * (1 / rho[i] - 1), 1), gamma2((1 - p[i]) * (1 / rho[i] - 1), 1);
            binomial_distribution<int> bino(n[i], 1 / (1 + gamma2(rng) / gamma1(rng)));
            r[i] = bino(rng);
        }
    }
    
    
    return;
}


