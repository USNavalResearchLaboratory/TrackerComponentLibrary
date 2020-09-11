/**ASSIGN3DLB Obtain a lower bound for the cost value of an axial 3D
 *           assignment problem. Such problems are NP-hard and thus cannot
 *           be solved in polynomial time. The optimization problem being
 *           bounded is
 *           minimize
 *           sum_{i=1}^{n1}sum_{j=1}^{n2}sum_{k=1}^{n3}C_{i,j,k}*rho_{i,j,k}
 *           subject to
 *           sum_{i=1}^{n1}sum_{j=1}^{n2}rho_{i,j,k}<=1 for all k
 *           sum_{i=1}^{n1}sum_{k=1}^{n3}rho_{i,j,k}<=1 for all j
 *           sum_{j=1}^{n2}sum_{k=1}^{n3}rho_{i,j,k} =1 for all i
 *           rho_{i,j,k} = 0 or 1
 *           assuming that n1<=n2<=n3, and C is an n1Xn2Xn3 cost matrix.
 *           Additional comments and examples are given in the Matlab
 *           implementation of the function.
 *
 *INPUTS: C An n1Xn2Xn3 cost hypermatrix. The costs are real numbers >-Inf.
 *          n1<=n2<=n3; If an empty matrix is passed, then lowerBound=0 is
 *          returned.
 *  method An optional parameter selecting the bound algorithm to use.
 *         Possible values are:
 *         0 Use the projection method followed by the Hungarian algorithm
 *           as in [1] (implemented using the Jonker-Volgenant algorithm in
 *           assign2D in place of the Hungarian algorithm). This algorithm
 *           requires that n1=n2=n3, so if that is not the case, the cost
 *           matrix is implicitly augmented so the lower bound can still be
 *           used. The implicit augmentation is discussed in the comments
 *           to the code.
 *         1 Use the simple (first) method of summing the minimum values
 *           across each first index of C as in [2]. Note that permuting
 *           the indices of C can change the tightness of the bound and the
 *           speed of this approach.
 *         2 (The default if omitted or an empty matrix is passed) Use the
 *           lower bound from the dual cost of the Lagrangian relaxation
 *           down to 2D assignment, computed setting the dual variables all
 *           to zero.
 *
 *OUTPUTS: lowerBound A lower bound on the value of the axial 3D assignment
 *                    optimization problem.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * lowerBound=assign3DLB(C,method)
 *
 *REFERENCES:
 *[1] B.-J. Kim, W. L. Hightower, P. M. Hahn, Y.-R. Zhu, and L. Sun, "Lower
 *    bounds for the axial three-index assignment problem," European
 *    Journal of Operational Research, vol. 202, no. 3, pp. 654-668, May
 *    2010.
 *[2] W. P. Pierskalla, "The tri-substitution method for the three-
 *   dimensional assignment problem," Canadian Operational Research Society
 *   Journal, vol. 5, no. 2, pp. 71-81, Jul. 1967.
 *
 *March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"

//Prototypes for the actual bound functions.
#include "assignAlgs3D.h"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    int method=2;
    size_t nVals[3];
    const double *C;//Initialized later from prhs[0].
    double lowerBound;
    
    if(nrhs>2||nrhs<1){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }

    //Empty matrix special case.
    if(mxIsEmpty(prhs[0])) {
        plhs[0]=mxCreateDoubleScalar(0);
        return;
    }
    
    checkRealDoubleHypermatrix(prhs[0]);
    {
        const mwSize numDims=mxGetNumberOfDimensions(prhs[0]);
        const mwSize *theDims=mxGetDimensions(prhs[0]);
        
        if(numDims>3) {
            mexErrMsgTxt("C has too many dimensions.");
            return;
        }
        
       nVals[0]=theDims[0];
       nVals[1]=theDims[1];
        
        //numDims will never be <2 due to how Matlab parameterizes things.
        if(numDims==2) {
            nVals[2]=1;
        } else {
            nVals[2]=theDims[2];
        }
    }

    if(!(nVals[0]<=nVals[1]&&nVals[1]<=nVals[2])) {
        mexErrMsgTxt("It is required that n1>=n2>=n3");
        return;
    }

    C=mxGetDoubles(prhs[0]);
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
         method=getIntFromMatlab(prhs[1]);
    }
    
    switch(method) {
        case 0://The projection-Hungarian algorithm of [1].
        {
            const size_t bufferSize=assign3DLBHungarianBufferSize(nVals);
            void * tempSpace=mxMalloc(bufferSize);
            
            lowerBound=assign3DLBHungarian(nVals,C, tempSpace);
            mxFree(tempSpace);
        }
            break;
        case 1://The very simple (first) method of [2].
            lowerBound=assign3DLBCPierskalla(nVals,C);
            break;
        case 2://The dual cost with zero dual variables in [3].
        {
            const size_t bufferSize=assign3DLBDual0BufferSize(nVals);
            void * tempSpace=mxMalloc(bufferSize);
            
            lowerBound=assign3DLBDual0(nVals,C, tempSpace);
            mxFree(tempSpace);
        }
            break;
        default:
            mexErrMsgTxt("Invalid method specified.");
            return;
    } 
    
    plhs[0]=mxCreateDoubleScalar(lowerBound);
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
