/**FACTORIALLOOKUP Evaluate a factorial using a lookup table. An
*  implementation in C++. See the comments to the Matlab implementation for
*  details on the inputs.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* val=factorialLookup(n);
*
*October 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
/*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"
#include "mathFuncs.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs<1){
        mexErrMsgTxt("Not enough inputs.");
        return;
    }
    
    if(nrhs>1) {
        mexErrMsgTxt("Too many inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Too many outputs.");
    }

    checkRealDoubleArray(prhs[0]);
    const size_t M=mxGetM(prhs[0]);
    const size_t N=mxGetN(prhs[0]);
    const size_t numEls=M*N;

    if(numEls==0) {
        plhs[0]=mxCreateDoubleMatrix(M,N,mxREAL);
        return;
    }

    const double *x=mxGetDoubles(prhs[0]);
    mxArray *retMat=mxCreateDoubleMatrix(M,N,mxREAL);
    double * retVals=mxGetDoubles(retMat);

    for(size_t i=0;i<numEls;i++) {
        retVals[i]=factorialLookupCPP(x[i]);
    }

    plhs[0]=retMat;
}

/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
