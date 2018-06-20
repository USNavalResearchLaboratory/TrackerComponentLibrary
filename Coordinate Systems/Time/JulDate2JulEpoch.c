/**JULDATE2JULEPOCH Convert two-part Julian dates given in a uniform
 *                  timescale, such as terrestrial time (TT) date to a
 *                  Julian epoch in the same scale. The timescale
 *                  should be uniform (i.e. not UTC).
 *
 *INPUTS: Jul1, Jul2  Matrices of two parts of a Julian date given in the 
 *                    same timescale as the epoch The units of the date are
 *                    days. The full date is the sum of both terms. The
 *                    entries in the matrices correspond to different
 *                    dates to be converted.
 *
 *OUTPUTS: epochJ A matrix of Julian epochs as fractional years in the same
 *                timescale as the input.
 *
 *A Julian epoch is a factional year number denominated in terms of a
 *year of exactly 365.25 days in TT. For example 2000.5 is a Julian epoch
 *in the middle of the year 2000.
 *
 *This is a wrapper for the function iauEpj in the International
 *Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
 *
 *The algorithm can be compiled for use in Matlab  using the
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *epochJ=JulDate2JulEpoch(Jul1,Jul2);
 *
 *March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

 /*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *Jul1,*Jul2,*epochJ;
    size_t numRow, numCol, numElements,i;
    mxArray *JulianEpochRetMATLAB;
    
    if(nrhs!=2){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>1){
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    numRow=mxGetM(prhs[0]);
    numCol=mxGetN(prhs[0]);
    numElements=numRow*numCol;
    
    if(numElements==0||numRow!=mxGetM(prhs[1])||numCol!=mxGetN(prhs[1])) {
        mexErrMsgTxt("The dimensionalities of the inputs are incorrect.");
        return;
    }
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    Jul1=(double*)mxGetData(prhs[0]);
    Jul2=(double*)mxGetData(prhs[1]);
    
    JulianEpochRetMATLAB=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    epochJ=(double*)mxGetData(JulianEpochRetMATLAB);

    for(i=0;i<numElements;i++) {
    /*Call the function in the SOFA library.*/
        epochJ[i]=iauEpj(Jul1[i], Jul2[i]);
    }
    
    plhs[0]=JulianEpochRetMATLAB;
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
