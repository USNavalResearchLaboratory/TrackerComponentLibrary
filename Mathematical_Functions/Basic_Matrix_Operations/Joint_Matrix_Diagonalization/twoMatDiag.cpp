/**TWOMATDIAG Given two real or complex Hermitian matrices, the first of
 *     which must be positive definite, find a matrix W that diagonalizes
 *     both of them. Specifically, W*C1*W'=I and W*C2*W'=a diagonal matrix.
 *     See the comments in the Matlab implementation of this file for more
 *     information on the inputs to the function and how it works.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * W=twoMatDiag(C1,C2,algorithm)
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
    mxArray *ComplexifiedArray=nullptr;
    std::complex <double> *C1Complex, *C2Complex;
    double *C1Real, *C2Real;
    bool isComplex;
    
    if(nrhs>3||nrhs<2){
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
    if(N!=mxGetN(prhs[0])||N!=mxGetM(prhs[1])||N!=mxGetN(prhs[1])) {
        mexErrMsgTxt("The matrix dimensions are inconsistent.");
        return;
    }
    
    if(nrhs>2) {
        algorithm=getIntFromMatlab(prhs[2]);
    }
    
    if(mxIsComplex(prhs[0])&&mxIsComplex(prhs[1])) {
        C1Complex=reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[0]));
        C2Complex=reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[1]));
        isComplex=true;
    } else if(mxIsComplex(prhs[0])&&!mxIsComplex(prhs[1])) {
        C1Complex=reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[0]));
        ComplexifiedArray=mxDuplicateArray(prhs[1]);
        mxMakeArrayComplex(ComplexifiedArray);
        C2Complex=reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(ComplexifiedArray));
        
        isComplex=true;
    } else if(!mxIsComplex(prhs[0])&&mxIsComplex(prhs[1])) {
        C2Complex=reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[1]));
        ComplexifiedArray=mxDuplicateArray(prhs[0]);
        mxMakeArrayComplex(ComplexifiedArray);
        C1Complex=reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(ComplexifiedArray));

        isComplex=true;
    } else {   
        //Both arrays are real.
        C1Real=mxGetDoubles(prhs[0]);
        C2Real=mxGetDoubles(prhs[1]);
        isComplex=false;
    }

    mxArray *WMatlab;
    if(isComplex) {
        const Eigen::Map<Eigen::MatrixXcd> C1Eigen(C1Complex,N,N);
        const Eigen::Map<Eigen::MatrixXcd> C2Eigen(C2Complex,N,N);
        
        Eigen::MatrixXcd WEigen(N,N);
        switch(algorithm) {
            case 0:
                WEigen=twoMatDiagSVD<Eigen::MatrixXcd>(C1Eigen,C2Eigen);
                break;
            case 1:
                WEigen=twoMatDiagEig<Eigen::MatrixXcd>(C1Eigen,C2Eigen);
                break;
            default:
                mexErrMsgTxt("Unknown algorithm specified.");
                return;
        }
        
        WMatlab=mxCreateDoubleMatrix(N,N,mxCOMPLEX);
        std::complex<double>* W=reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(WMatlab));
        for(size_t k=0;k<N*N;k++) {
            W[k]=WEigen(k);
        }
        
    }  else {
        const Eigen::Map<Eigen::MatrixXd> C1Eigen(C1Real,N,N);
        const Eigen::Map<Eigen::MatrixXd> C2Eigen(C2Real,N,N);
        
        Eigen::MatrixXd WEigen(N,N);
        switch(algorithm) {
            case 0:
                WEigen=twoMatDiagSVD<Eigen::MatrixXd>(C1Eigen,C2Eigen);
                break;
            case 1:
                WEigen=twoMatDiagEig<Eigen::MatrixXd>(C1Eigen,C2Eigen);
                break;
            default:
                mexErrMsgTxt("Unknown algorithm specified.");
                return;
        }

        WMatlab=mxCreateDoubleMatrix(N,N,mxREAL);
        double *W=mxGetDoubles(WMatlab);
        for(size_t k=0;k<N*N;k++) {
            W[k]=WEigen(k);
        }
    }

    //The return value.
    plhs[0]=WMatlab;
    
    if(ComplexifiedArray!=nullptr) {
        mxDestroyArray(ComplexifiedArray);
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
 