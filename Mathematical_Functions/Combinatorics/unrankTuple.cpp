/**UNRANKTUPLE Obtain the tuple corresponding to its order in the sequence
 *            of tuples given a certain set of maximum values for each
 *            digit. These are generated with the first digit being the
 *            LEAST significant. See the comments to the Matlab version for
 *            more information.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *tuple=unrankTuple(rank,maxVals,areProdVals)
 *
 *November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"

#include "combinatorialFuns.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    bool areProdVals=false;

    if(nrhs<2||nrhs>3) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    if(nlhs>1) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;  
    }

    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        areProdVals=getBoolFromMatlab(prhs[2]);
    }

    //If an empty matrix is passed, return an empty matrix.
    if(mxIsEmpty(prhs[1])) {
        plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
        return;
    }

    //The default datatype in Matlab is double, so we assume that is what
    //the inputs and outputs will be. Things are slower if we allocate
    //space and copy these arrays rather than if we typecast stuff. 
    size_t numTuples=mxGetNumberOfElements(prhs[0]);
    checkRealDoubleArray(prhs[0]);
    const double *tuples=mxGetDoubles(prhs[0]);
    size_t numDimEls=mxGetNumberOfElements(prhs[1]);
    checkRealDoubleArray(prhs[1]);
    const double *dims=mxGetDoubles(prhs[1]);
    
    for(size_t i=0;i<numDimEls;i++) {
        if(dims[i]<1) {
            mexErrMsgTxt("The elements of the maxVals input must be >=1.");
            return; 
        }
    }
    
    mxArray *newTupleMat;
    size_t numDim;
    bool isPastEnd;
    if(areProdVals==true) {
        numDim=numDimEls-1;
        const double *maxVals=dims;
        newTupleMat=mxCreateDoubleMatrix(numDim,numTuples,mxREAL);
        double *newTuples=mxGetDoubles(newTupleMat);
        
        isPastEnd=unrankFloatTupleCPP(numDim, numTuples, tuples, maxVals, newTuples);
    } else {
        numDim=numDimEls;
        double *tempSpace=new double[numDim+1];
        newTupleMat=mxCreateDoubleMatrix(numDim,numTuples,mxREAL);
        double *newTuples=mxGetDoubles(newTupleMat);
        
        isPastEnd=unrankFloatTupleDimsCPP(numDim, numTuples, tuples, dims, newTuples, tempSpace);
           
        delete[] tempSpace;
    }
    
    //Return an empty matrix if an index is beyond the end of the array.
    if(isPastEnd) {
        mxDestroyArray(newTupleMat);
        plhs[0]=mxCreateDoubleMatrix(numDim,0,mxREAL);
        return;
    }
    
    plhs[0]=newTupleMat;
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
