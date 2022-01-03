/**BACKSUBSTITUTION Given an upper-triangular matrix, use back substitution
*       to solve the system U*x=b for x.  See the Matlab implementation for
*       more details. This is only implemented for real matrices.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* b=backSubstitution(U,b,algorithm);
*
*June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"

/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "basicMatOpsEigen.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    int algorithm=0;
    
    if(nrhs<2||nrhs>3){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>1) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }
    
    checkDoubleArray(prhs[0]);
    checkDoubleArray(prhs[1]);
    const size_t N=mxGetM(prhs[0]);
    
    if(mxGetM(prhs[0])!=N) {
        mexErrMsgTxt("The L matrix must be square.");
        return;
    }
    
    if(mxGetM(prhs[1])!=N||mxGetN(prhs[1])!=1) {
        mexErrMsgTxt("The dimensions b are incorrect.");
        return;
    }

    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        algorithm=getIntFromMatlab(prhs[2]);
    }
    
    double *L=mxGetDoubles(prhs[0]);
    const double *b=mxGetDoubles(prhs[1]);
    const Eigen::Map<Eigen::MatrixXd> UEigen(L,N,N);

    Eigen::MatrixXd bEigen(N,1);
    for(size_t k=0;k<N;k++) {
        bEigen(k)=b[k];
    }

    switch(algorithm) {
        case 0:
            backSubstitutionByRow<Eigen::MatrixXd>(UEigen,bEigen);
            break;
        case 1:
            backSubstitutionByCol<Eigen::MatrixXd>(UEigen,bEigen);
            break;
        default:
            mexErrMsgTxt("Invalid algorithm specified.");
            return;
    }

    //For the return values.
    mxArray *retVec=mxCreateDoubleMatrix(N,1,mxREAL);
    double *retBuff=mxGetDoubles(retVec);
    
    for(size_t k=0;k<N;k++) {
        retBuff[k]=bEigen(k);
    }
    
    plhs[0]=retVec;
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
 