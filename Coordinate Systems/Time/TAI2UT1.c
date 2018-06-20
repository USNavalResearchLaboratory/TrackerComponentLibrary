/**TAI2UT1 Convert from international atomic time (TAI) to UT1, which is a
 *         nonuniform timescale based on the rotation rate of the Earth.
 *
 *INPUTS: Jul1,Jul2 Two parts of a pseudo-Julian date given in TAI. The
 *                  units of the date are days. The full date is the sum of
 *                  both terms. The date is broken into two parts to
 *                  provide more bits of precision. It does not matter how
 *                  the date is split.
 *           deltaT An optional parameter specifying the offset between TAI
 *                  and UT1 in seconds. If this parameter is omitted, then
 *                  a value derived from function getEOP will be used.
 *
 *OUTPUTS: Jul1,Jul2 Two parts of a Julian date in UT1.
 *
 *This is a wrapper for the function iauTaiut1 in the International
 *Astronomical Union's Standards of Fundamental Astronomy library.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[Jul1,Jul2]=TAI2UT1(Jul1,Jul2);
 *
 *Many temporal coordinate systems standards are compared in [1].
 *
 *REFERENCES:
 *[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
 *    Temporal Coordinate Systems for Target Tracking," Formal Report,
 *    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
 *    173 pages.
 *
 *March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double Jul1, Jul2,deltaT;
    
    int retVal;
    
    if(nrhs<2||nrhs>3){
        mexErrMsgTxt("Wrong number of inputs.");
        return;
    }

    if(nlhs>2){
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    Jul1=getDoubleFromMatlab(prhs[0]);
    Jul2=getDoubleFromMatlab(prhs[1]);
    
    if(nrhs>2) {
        deltaT=getDoubleFromMatlab(prhs[2]);
    } else {
        mxArray *retVals[4];
        mxArray *JulUTCMATLAB[2];
        double JulUTC[2];
        
        //Get the time in UTC to look up the parameters.
        retVal=iauTaiutc(Jul1, Jul2, &JulUTC[0], &JulUTC[1]);
        switch(retVal){
            case 1:
                mexWarnMsgTxt("Dubious Date entered.");
                break;
            case -1:
                mexErrMsgTxt("Unacceptable date entered");
                break;
            default:
                break;
        }
        
        JulUTCMATLAB[0]=doubleMat2Matlab(&JulUTC[0],1,1);
        JulUTCMATLAB[1]=doubleMat2Matlab(&JulUTC[1],1,1);
        
        //Get the Earth orientation parameters for the given date.
        mexCallMATLAB(4,retVals,2,JulUTCMATLAB,"getEOP");
        //The 32.184 is the offset between TT and TAI.
        deltaT=getDoubleFromMatlab(retVals[3])-32.184;
        //Free the returned arrays.
        mxDestroyArray(retVals[0]);
        mxDestroyArray(retVals[1]);
        mxDestroyArray(retVals[2]);
        mxDestroyArray(retVals[3]);
        mxDestroyArray(JulUTCMATLAB[0]);
        mxDestroyArray(JulUTCMATLAB[1]);
    }
 
    //Perform the conversion.
    retVal=iauTaiut1(Jul1, Jul2, -deltaT, &Jul1, &Jul2 );
    switch(retVal) {
        case -1:
            mexErrMsgTxt("Unacceptable date provided.");
            return;
        case 1:
            mexWarnMsgTxt("Dubious year provided\n");
            break;
        default:
            break;
    }
    
    plhs[0]=doubleMat2Matlab(&Jul1,1, 1);
    if(nlhs>1) {
        plhs[1]=doubleMat2Matlab(&Jul2,1, 1);
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
