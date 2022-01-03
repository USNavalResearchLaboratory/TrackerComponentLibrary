/**PESUDOINVERSE Compute the matrix pseudoinverse using a singular value
*               decomposition. See the Matlab implementation for more
*               details.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* pinvX=pseudoInverse(X,algorithm)
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
    
    if(nrhs<1||nrhs>2){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>1) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }
    
    if(nrhs>1) {
        algorithm=getIntFromMatlab(prhs[1]);
    }
 
    checkDoubleArray(prhs[0]);
    
    const size_t M=mxGetM(prhs[0]);
    const size_t N=mxGetN(prhs[0]);
    
    if(mxIsComplex(prhs[0])) {
        std::complex<double> *X=reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[0]));
        const Eigen::Map<Eigen::MatrixXcd> XEigen(X,M,N);
        Eigen::MatrixXcd pinvX;

        pinvX=pseudoInverse<Eigen::MatrixXcd>(XEigen, algorithm);
        mxArray *retMat=mxCreateDoubleMatrix(N,M,mxCOMPLEX);
        std::complex<double> *pRet=reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(retMat));
        
        const size_t MN=M*N;
        for(size_t k=0;k<MN;k++) {
          pRet[k]=pinvX(k);  
        }
        
        plhs[0]=retMat;
    } else {
        double *X=mxGetDoubles(prhs[0]);
        const Eigen::Map<Eigen::MatrixXd> XEigen(X,M,N);
        Eigen::MatrixXd pinvX;
        
        pinvX=pseudoInverse<Eigen::MatrixXd>(XEigen, algorithm);
        mxArray *retMat=mxCreateDoubleMatrix(N,M,mxREAL);
        double *pRet=mxGetDoubles(retMat);
        
        const size_t MN=M*N;
        for(size_t k=0;k<MN;k++) {
          pRet[k]=pinvX(k);  
        }
        plhs[0]=retMat;
    }
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
 