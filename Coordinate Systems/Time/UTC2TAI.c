/**UTC2TAI Convert from universal coordinated time (UTC) given as a two-
 *         part pseudo-Julian date to international atomic time (TAI),
 *         represented as a two-part pseudo-Julian date.
 *
 *INPUTS: Jul1, Jul2 Matrices of two parts of a pseudo-Julian date given in
 *                   UTC. The units of the date are days. The full date is
 *                   the sum of both terms. The date is broken into two
 *                   parts to provide more bits of precision. It does not
 *                   matter how the date is split. Corresponding elements
 *                  in each matrix are times that are converted.
 *
 *OUTPUTS: Jul1, Jul2 The time as a Julian date in TAI with the same
 *                     dimensionalities as the input sets of dates.
 *   
 *This is a mex wrapper for the function iauUtctai in the International
 *Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
 *The UTC date is only pseudo-Julian, because there is not a fixed number
 *of seconds in a Julian day. The convention used in the IAU standard is
 *that the Julian day matches the UTC day regardless of whether the UTC day
 *is 86399, 86400 or 86401 SI seconds (depending on the presence of leap
 *seconds).
 *
 *UTC began at 1960 January 1.0 (JD 2436934.5) and this function should not
 *be called with an earlier date.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[Jul1,Jul2]=UTC2TAI(Jul1,Jul2);
 *
 *Many temporal coordinate systems standards are compared in [1].
 *
 *REFERENCES:
 *[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
 *    Temporal Coordinate Systems for Target Tracking," Formal Report,
 *    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
 *    173 pages.
 *
 * October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *Jul1,*Jul2, *Jul1Ret, *Jul2Ret;
    mxArray *Jul1RetMATLAB,*Jul2RetMATLAB;
    size_t numRow, numCol, numElements,i;
    
    if(nrhs!=2){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>2){
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
    
    //Allocate space for the return variables.
    Jul1RetMATLAB=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    Jul1Ret=(double*)mxGetData(Jul1RetMATLAB);
    Jul2RetMATLAB=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    Jul2Ret=(double*)mxGetData(Jul2RetMATLAB);
    
    for(i=0;i<numElements;i++) {
        int retVal;
    /*Call the function in the SOFA library.*/
        retVal=iauUtctai(Jul1[i], Jul2[i], Jul1Ret+i, Jul2Ret+i);

        switch(retVal){
            case 1:
                mexWarnMsgTxt("Dubious Date entered.");
                break;
            case -1:
                mxDestroyArray(Jul1RetMATLAB);
                mxDestroyArray(Jul2RetMATLAB);
                mexErrMsgTxt("Unacceptable date entered");
                break;
            default:
                break;
        }
    }
    
    plhs[0]=Jul1RetMATLAB;
    if(nlhs>1) {
        plhs[1]=Jul2RetMATLAB;
    } else {
        mxDestroyArray(Jul2RetMATLAB);
    }
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
