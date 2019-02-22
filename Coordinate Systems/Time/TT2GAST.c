/**TT2GAST Convert from terrestrial time (TT) to Greenwhich apparent 
 *         sidereal time (GMST), which is a measure of the rotational
 *         direction of the Earth.
 *
 *INPUTS: Jul1,Jul2 Two parts of a pseudo-Julian date given in TT. The
 *                  units of the date are days. The full date is the sum of
 *                  both terms. The date is broken into two parts to
 *                  provide more bits of precision. It does not matter how
 *                  the date is split.
 *          version An optional integer specifying the theory to use for
 *                  GAST. The theory chosen should be consistent with other
 *                  values used in astronomical routines. Possible values
 *                  are
 *                  1994 Compute GAST ion accordance with the International
 *                     Astronomical Union's (IAU's) 1994 model.
 *                  2000 Compute GAST in line with IAU 2000 resolutions
 *                     related to precession and nutation.
 *                  2006 (The default if omitted) Compute GAST in line with
 *                     IAU 2006 resolutions related to precession and
 *                     nutation.
 *           deltaT An optional parameter specifying the offset between TT
 *                  and UT1 in seconds. If this parameter is omitted, then
 *                  the value of the function getEOP will be used.
 *
 *OUTPUTS: GAST The Greenwhich apparent sideral time in radians. 
 *
 *GAST is defined in Section 5.5.7 of [1].
 *
 *This is a wrapper for the functions iauTtut1 and iauGst94, iauGst00a,
 *and iauGst06a in the International Astronomical Union's Standards of
 *Fundamental Astronomy library.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *GAST=TT2GAST(Jul1,Jul2);
 *or
 *GAAST=TT2GAST(Jul1,Jul2,version);
 *or
 *GAST=TT2GAST(Jul1,Jul2,version,deltaT);
 *
 *REFERENCES:
 *[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
 *    Rotation and Reference Systems Service Std. 36, 2010.
 *
 *April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double TT1, TT2, UT11,UT12, deltaT,GAST;
    int algToUse=2006;
    int retVal;
    
    if(nrhs<2||nrhs>4){
        mexErrMsgTxt("Wrong number of inputs.");
        return;
    }

    if(nlhs>1){
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    TT1=getDoubleFromMatlab(prhs[0]);
    TT2=getDoubleFromMatlab(prhs[1]);
    
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        algToUse=getIntFromMatlab(prhs[2]);
    }
    
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        deltaT=getDoubleFromMatlab(prhs[3]);
    } else {
        mxArray *retVals[4];
        mxArray *JulUTCMATLAB[2];
        double JulUTC[2];
        
        //Get the time in UTC to look up the parameters by going to TAI and
        //then UTC.
        retVal=iauTttai(TT1, TT2, &JulUTC[0], &JulUTC[1]);
        if(retVal!=0) {
            mexErrMsgTxt("An error occurred computing TAI.");
        }
        
        retVal=iauTaiutc(JulUTC[0], JulUTC[1], &JulUTC[0], &JulUTC[1]);
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
        //This is TT-UT1
        deltaT=getDoubleFromMatlab(retVals[3]);
        //Free the returned arrays.
        mxDestroyArray(retVals[0]);
        mxDestroyArray(retVals[1]);
        mxDestroyArray(retVals[2]);
        mxDestroyArray(retVals[3]);
        mxDestroyArray(JulUTCMATLAB[0]);
        mxDestroyArray(JulUTCMATLAB[1]);
    }
     
    //Get UT1
    retVal=iauTtut1(TT1, TT2, deltaT, &UT11, &UT12);
    if(retVal!=0) {
        mexErrMsgTxt("An error occurred computing UT1.");
    }
    
    //Get Greenwhich mean sidereal time in radians using the chosen
    //algorithm
    switch(algToUse) {
        case 1994:
            GAST=iauGst94(UT11, UT12);
            break;
        case 2000:
            GAST=iauGst00a(UT11, UT12, TT1, TT2);
            break;
        case 2006:
            GAST=iauGst06a(UT11, UT12, TT1, TT2);
            break;
        default:
            mexErrMsgTxt("An invalid algorithm version was given.");
    }

    plhs[0]=doubleMat2Matlab(&GAST,1, 1);
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
