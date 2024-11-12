/**BINOMIAL A C++ implementation of a function to evaluate the binomial
*          coefficient n choose k, the total number of ways of choosing k
*          items from a total selection of n where the ordering of the
*          selection does not matter. See the Matlab implementation for
*          more comments.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* val=binomial(n,k);
*
*September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
*Interface translated into C++:
*September 2024 Kathleen Appenzeller, WR Systems, Fairfax VA */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release. */

#include "MexValidation.h"
/*This header is required by Matlab*/
#include "mex.h"
#include "combinatorialFuns.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {

    if(nrhs<1){
        mexErrMsgTxt("Not enough inputs.");
    }
    
    if(nrhs>2){
        mexErrMsgTxt("Too many inputs.");
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Too many outputs.");
    }

    if (mxIsEmpty(prhs[0])){
        mexErrMsgTxt("Missing input argument n");
    }
    const size_t n = getSizeTFromMatlab(prhs[0]);

    if(n>9007199254740992) {//That number is 2^53
        mexErrMsgTxt("n must be less than or equal to 2^53.");
    }
    if (mxIsEmpty(prhs[1])){
        mexErrMsgTxt("Missing input argument k");
    }

    const size_t k = getSizeTFromMatlab(prhs[1]);
    if(k>9007199254740992) {//That number is 2^53.
        mexErrMsgTxt("k must be less than or equal to 2^53.");
    }

    const double val = binomialCPP(n, k);

    mxArray * retIdx = mxCreateDoubleScalar(val);
    plhs[0]=retIdx;
}