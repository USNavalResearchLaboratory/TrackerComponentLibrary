/**ASSIGN3DBB Solve the data fusion axial 3D assignment problem using a
 *        branch-and-bound algorithm. Such problems are NP-hard and thus
 *        cannot be solved in polynomial time, so this function can only
 *        solve problems of a moderate size. See the comments in the Matlab
 *        implementation of this file for more information on the inputs to
 *        the function and how it works.
 *
 * It is assumed that all matrices are sufficiently small that it does not
 * matter whether lengths are held in size_t or ptrdiff_t variables.
 *
 * The algorithm can be compiled for use in Matlab  using the 
 * CompileCLibraries function.
 *
 * The algorithm is run in Matlab using the command format
 * [minCostTuples,gain]=assign3DBB(C,maximize,boundType,initMethod,maxIter,epsVal)
 *
 *October 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"

/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"

//Prototypes for the actual assignment functions.
#include "assignAlgs3DCPP.hpp"
#include <limits>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    bool maximize=false;
    int boundType=2;
    int initMethod=0;
    size_t maxIter=10;
    double epsVal=std::numeric_limits<double>::epsilon();
    size_t nDims[3];

    if(nrhs>6||nrhs<1){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>2) {
        mexErrMsgTxt("Invalid number of outputs.");
        return;
    }

    ////////
    //Check the matrix C
    ///////
    if(mxIsComplex(prhs[0])) {
        mexErrMsgTxt("C must be real.");
        return;
    }
    
    if(mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) {
        mexErrMsgTxt("C must be of type double.");
        return;
    }

    if(mxIsEmpty(prhs[0])) {//The empty matrix special case.
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);//tuples=[];
        
        if(nlhs>1) {
            plhs[1]=mxCreateDoubleScalar(0);//gain=0
        }
        return;
    }
    
    const double *const COrig=mxGetDoubles(prhs[0]);    
    if(mxGetNumberOfDimensions(prhs[0])!=3&&!mxIsScalar(prhs[0])) {
        mexErrMsgTxt("C must be 3D.");
        return;
    }
    
    //The algorithm was implemented assuming that
    //nDims[0]<=nDims[1]<=nDims[2] for the
    //dimensionality of C. If a matrix having a different arrangement of
    //indices is passed, one could permute the indices to make the
    //assumptions below hold.
    {
        const size_t *nDimsMatlab=mxGetDimensions(prhs[0]);

        //Deal with trailing singleton dimensions.
        switch(mxGetNumberOfDimensions(prhs[0])) {
            case 3:
                nDims[0]=nDimsMatlab[0];
                nDims[1]=nDimsMatlab[1];
                nDims[2]=nDimsMatlab[2];
                break;
            case 2:
                nDims[0]=nDimsMatlab[0];
                nDims[1]=nDimsMatlab[1];
                nDims[2]=1;
                break;
            case 1:
                nDims[0]=nDimsMatlab[0];
                nDims[1]=1;
                nDims[2]=1;
                break;
            default:
                mexErrMsgTxt("C has an incorrect number of dimensions.");
        }
    }
    
    if(!(nDims[0]<=nDims[1]&&nDims[1]<=nDims[2])) {
        mexErrMsgTxt("It is required that size(C,1)<=size(C,2)<=size(C,3)");
        return;
    }

    ////////
    //Check the other inputs
    ///////
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        maximize=getBoolFromMatlab(prhs[1]);
    }
    
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        boundType=getIntFromMatlab(prhs[2]);
        
        if(boundType<0||boundType>2) {
            mexErrMsgTxt("Invalid boundType specified.");
            return;
        }
    }
    
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        initMethod=getIntFromMatlab(prhs[3]);
        
        if(initMethod<0||initMethod>1) {
            mexErrMsgTxt("Invalid initMethod specified.");
            return;
        }
    }
    
    if(nrhs>4&&!mxIsEmpty(prhs[4])) {
        maxIter=getSizeTFromMatlab(prhs[4]);
        
        if(maxIter<1) {
            mexErrMsgTxt("maxIter must be >=1.");
            return;
        }
    }
    
    if(nrhs>5&&!mxIsEmpty(prhs[5])) {
        epsVal=getDoubleFromMatlab(prhs[5]);
        
        if(maxIter<0) {
            mexErrMsgTxt("epsVal must be >=0.");
            return;
        }
    }
    
    //Allocate space for a return value.
    ptrdiff_t *minCostTuples=new ptrdiff_t[3*nDims[0]];
    
	const size_t bufferSize=assign3DBBBufferSize(nDims,boundType,initMethod);
    uint8_t *tempBuffer=new uint8_t[bufferSize];
    
    const double gain=assign3DBBCPP(minCostTuples,
                                            COrig,
                                            nDims,
                                         maximize,
                                        boundType,
                                       initMethod,
                                          maxIter,
                                           epsVal,
              reinterpret_cast<void*>(tempBuffer));
    
    delete[] tempBuffer;

    //If no feasible solution was found.
    if(!std::isfinite(gain)) {
        //Return an empty matrix.
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
    } else {
        //Convert the indexation to Matlab indexation.
        for(size_t k=0;k<3*nDims[0];k++) {
            minCostTuples[k]++;
        }
        
        //Set the output.
        plhs[0]=ptrDiffTMat2MatlabDoubles(minCostTuples,3,nDims[0]);
    }
    
    delete[] minCostTuples;
    
    if(nlhs>1) {
        plhs[1]=mxCreateDoubleScalar(gain);
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
*OF RECIPIENT IN THE USE OF THE SOFTWARE.
*/
