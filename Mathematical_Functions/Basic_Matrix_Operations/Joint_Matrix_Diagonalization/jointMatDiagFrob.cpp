/**JOINTMATDIAGFROB Given K real NXN symmetric matrices in A, try to find a
*               common diagonalization matrix V such that V*A(:,:,k)*V' is
*               diagonal for all k. V is generally not orthogonal. When
*               K>2, it is possible that no common matrix exists, in which
*               case, the sum of the squared Frobenius norms of the
*               matrices E_k for k=1:K is minimized. E_k is the matrix 
*               V*A(:,:,k)*V' with the main diagonal elements set to zero.
*               Alternatively, positive weights can be provided in w to
*               allow for the Frobenius norms to be weighred going into the
*               cost function. See the comments to the Matlab
*               implementation for more information.
*   
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* [V,C,FCost,exitCode]=jointMatDiagFrob(A,w,maxIter,RelTol,AbsTol)
*
*June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"

/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "basicMatOpsEigen.hpp"
#include <limits>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t maxIter=1000;
    double RelTol=std::numeric_limits<double>::epsilon();
    double AbsTol=std::numeric_limits<double>::epsilon();
    double *A;
    double *w=nullptr;
    size_t N, K=0;
    mxArray *VMat, *CMat, *FCostMat, *exitCodeMat;
    
    if(nrhs>5||nrhs<1){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>4) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }
    
    if(nrhs>4&&!mxIsEmpty(prhs[4])) {
        AbsTol=getDoubleFromMatlab(prhs[4]);
    }
    
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        RelTol=getDoubleFromMatlab(prhs[3]);
    }
    
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        maxIter=getSizeTFromMatlab(prhs[2]);
    }
    
    checkRealDoubleHypermatrix(prhs[0]);
    N=mxGetM(prhs[0]);
    
    {
        const size_t numDims=mxGetNumberOfDimensions(prhs[0]);
        
        if(numDims>3) {
            mexErrMsgTxt("The first input has too many dimensions.");
            return;
        } else if(numDims==3) {
            const size_t *dimsVals=mxGetDimensions(prhs[0]);
            
            if(dimsVals[1]!=N) {
                mexErrMsgTxt("The A matrices are not square.");
                return;
            }

            K=dimsVals[2];
        } else if(mxGetNumberOfDimensions(prhs[0])==2) {
            K=1;
        } else {
            mexErrMsgTxt("Inconsistent number of dimensions.");
            return;
        }
    }
    
    A=mxGetDoubles(prhs[0]);

    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        checkRealDoubleArray(prhs[1]);
        const size_t numWRows=mxGetM(prhs[1]);
        const size_t numWCols=mxGetN(prhs[1]);
        
        if(numWRows!=1&&numWCols!=1) {
            mexErrMsgTxt("w cannot be a matrix.");
            return;
        }
        
        if(numWRows!=K&&numWCols!=K) {
            mexErrMsgTxt("w, if provided, must be length K.");
            return;
        }
        
        w=mxGetDoubles(prhs[1]);
        
        //If all the w values are equal, then treat it as there being no w.
        bool allEqual=true;
        for (size_t i=1;i<K;i++) {
            if(w[i]!=w[0]) {
                allEqual=false;
                break;
            }
        }
        
        if(allEqual) {
            w=nullptr;
        }
    }
    
    //Allocate the return matrices.
    VMat=mxCreateDoubleMatrix(N,N,mxREAL);
    {
        mwSize dims[3]={N,N,K};
        CMat=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    }
    FCostMat=mxCreateDoubleMatrix(1,1,mxREAL);
    exitCodeMat=mxCreateDoubleMatrix(1,1,mxREAL);
    
    double *V=mxGetDoubles(VMat);
    double *C=mxGetDoubles(CMat);
    double *FCost=mxGetDoubles(FCostMat);
    double *exitCode=mxGetDoubles(exitCodeMat);

    std::vector<Eigen::MatrixXd> AVec;
    //Make the vector able to hold all of the matrices.
    AVec.resize(K);
    //Put all of the slices of A into AVec.
    const size_t N2=N*N;
    for(size_t k=0;k<K;k++) {
        Eigen::Map<Eigen::MatrixXd> CurASlice(A+N2*k,N,N);
        AVec[k]=CurASlice;
    }

    //For the return values;
    Eigen::MatrixXd W(N,N);
    std::vector<Eigen::MatrixXd> CVec;
    
    if(w==nullptr) {
        //No weighting.
        *exitCode=static_cast<double>(JointMatDiagZieheAlg(AVec, maxIter, RelTol, AbsTol, W, CVec, *FCost));
    } else {
        //The matrices are weighted.
        Eigen::MatrixXd wEigen(N,1);
        //Normalize w.
        double wSum=0;
        for(size_t k=0;k<K;k++) {
           wSum+=w[k];
        }
        for(size_t k=0;k<K;k++) {
           wEigen(k)=w[k]/wSum;
        }
 
        *exitCode=static_cast<double>(JointMatDiagFrobVollAlg(AVec, wEigen, maxIter, RelTol, AbsTol,  W, CVec, *FCost));
    }

    //Copy the results into the return values.
    for(size_t i=0;i<N2;i++) {
        V[i]=W(i);
    }
    
    if(nlhs>1) {
        for(size_t k=0;k<K;k++) {
            const Eigen::MatrixXd curSlice=CVec[k];
            const size_t startIdx=N2*k;
            
            for(size_t i=0;i<N2;i++) {
                C[i+startIdx]=curSlice(i);
            }
        }
    }

    plhs[0]=VMat;
    if(nlhs>1) {
        plhs[1]=CMat;
        if(nlhs>2) {
            plhs[2]=FCostMat;
            if(nlhs>3) {
                plhs[3]=exitCodeMat;
            } else {
                mxDestroyArray(exitCodeMat);
            }
        } else {
            mxDestroyArray(FCostMat);
            mxDestroyArray(exitCodeMat);
        }
    } else {
        mxDestroyArray(CMat);
        mxDestroyArray(FCostMat);
        mxDestroyArray(exitCodeMat);
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
