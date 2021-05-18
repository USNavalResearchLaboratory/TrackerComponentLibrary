/*KRONSYM Take the symmetric Kronecker product of the matrices A and B. The
*         standard Kronecker product of two real, square matrices A and B,
*         one can write the following relation with respect to any real
*         square matrix S:
*         kron(A,B)*vec(S)==vec(B*S*A')
*         However, the symmetric Kronecker product is such that for any two
*         real square matrices A and B and a square SYMMETRIC matrix S, one
*         can write
*         kronSym(A,B)*vech(S,sqrt(2))==vech((1/2)*(B*S*A'+A*S*B'),sqrt(2))
*         Whereas is A and B are nXn, the standard kronecker product is
*         n^2Xn^2, the symmetric Kronecker product is
*         (n*(n+1)/2)X(n*(n+1)/2).
*
*INPUTS: A, B Two real nXn matrices. They need not be square.
*
*OUTPUTS: K The (n*(n+1)/2)X(n*(n+1)/2) symmetric Kronecker product of A
*           and B.
*
*Symmetric Kronecker products are discussed in the appendix of [1], where
*they play a role in the implementation of a semidefinite programming
*algorithm. More implementation details are given in the Matlab
*implementation of the algorithm.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*K=kronSym(A,B);
*
*REFERENCES:
*[1] F. Alizadeh, J.-P. A. Haeberly, and M. L. Overton, "Primal-dual
*    interior-point methods for semidefinite programming: Convergence
*    rates, stability and numerical results," SIAM Journal on Optimization,
*    vol. 8, no. 3, pp. 746-768, 1998.
*
*February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
//For sqrt
#include <math.h>
/*This is for input validation*/
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t n;
    size_t prodDim;
    double *A, *B, *K;
    mxArray *retMat;
    const double sqrt2=sqrt(2);
    const double dblSqrt2=2*sqrt2;
    
    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    if((mxIsEmpty(prhs[0])&&!mxIsEmpty(prhs[1]))||(!mxIsEmpty(prhs[0])&&mxIsEmpty(prhs[1]))) {
        mexErrMsgTxt("Invalid use of empty matrices as inputs.");
        return;
    }
    
    //If two empty matrices are passed, then just return an empty matrix.
    if(mxIsEmpty(prhs[0])&&mxIsEmpty(prhs[1])) {
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    
    if(mxGetNumberOfDimensions(prhs[0])>2||mxGetNumberOfDimensions(prhs[1])>2) {
        mexErrMsgTxt("This function does not work with hypermatrices.");
        return; 
    }
    
    n=mxGetN(prhs[0]);
    if(n!=mxGetM(prhs[0])||n!=mxGetN(prhs[1])||n!=mxGetM(prhs[1])) {
        mexErrMsgTxt("This function only works with pairs of squares matrices.");
        return;
    }
    
    prodDim=(n*(n+1))/2;
    
    //Allocate space for the return matrix
    retMat=mxCreateDoubleMatrix(prodDim,prodDim,mxREAL);
    A=mxGetDoubles(prhs[0]);
    B=mxGetDoubles(prhs[1]);
    K=mxGetDoubles(retMat);
    
    {
       size_t i1,i2,j1,j2;
        
        for(i1=1;i1<=n;i1++) {
            size_t col=i1+((i1-1)*(2*n-i1))/2;
            for(j1=1;j1<=n;j1++) {
                size_t KIdx=prodDim*(col-1)+j1+((j1-1)*(2*n-j1))/2-1;
                const size_t j1i1=j1+n*(i1-1)-1;
                size_t j2i1=j1i1;
                
                K[KIdx]=(A[j1i1]*B[j1i1]+B[j1i1]*A[j1i1])/2;
                for(j2=(j1+1);j2<=n;j2++) {
                    KIdx++;
                    j2i1++;
                    
                    K[KIdx]=(A[j1i1]*B[j2i1]+B[j1i1]*A[j2i1])/sqrt2;
                }
            }

            for(i2=(i1+1);i2<=n;i2++) {
                col++;
                for(j1=1;j1<=n;j1++) {
                    size_t KIdx=prodDim*(col-1)+j1+((j1-1)*(2*n-j1))/2-1;
                    const size_t j1i1=j1+n*(i1-1)-1;
                    const size_t j1i2=j1+n*(i2-1)-1;
                    size_t j2i1=j1i1;
                    size_t j2i2=j1i2;

                    K[KIdx]=(A[j1i1]*B[j1i2]+A[j1i2]*B[j1i1]+B[j1i1]*A[j1i2]+B[j1i2]*A[j1i1])/dblSqrt2;     
                    for(j2=(j1+1);j2<=n;j2++) {
                        KIdx++;
                        j2i1++;
                        j2i2++;
                        
                        K[KIdx]=(A[j1i1]*B[j2i2]+A[j1i2]*B[j2i1]+B[j1i1]*A[j2i2]+B[j1i2]*A[j2i1])/2;
                    }
                }
            }
        }
    }
    
    plhs[0]=retMat;
}

/*LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
