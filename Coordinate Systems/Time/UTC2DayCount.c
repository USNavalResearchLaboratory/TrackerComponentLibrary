/**UTC2DAYCOUNT A function for converting a universal coordinated time
 *              (UTC) given as a two-part pseudo-Julian to a date in terms
 *              of the year, day number and second in the day.
 *
 *INPUTS: Jul1, Jul2 Vectors or matrices of the two parts of Julian dates
 *                   given in UTC. The units of the date are days. The full
 *                   date is the sum of both terms. The date is broken into
 *                   two parts to provide more bits of precision. It does
 *                   not matter how the date is split.
 *
 *OUTPUTS: year A matrix of integer years (as Matlab doubles) in the
 *              Gregorian calendar under UTC corresponding to the input
 *              dates.
 *     dayCount A matrix integer days of the year (as Matlab doubles),
 *              accounting for leap years, corresponding to the input
 *              dates. Values in dayCount are integers >=1.
 *       second A matrix of the (possibly fractional) seconds of the day,
 *              accounting for leap seconds, corresponding to the input
 *              dates. Values in second are >=0.
 *
 *This makes use of the function iauJd2cal in the International
 *Astronomical Union's (IAU) Standard's of Fundamental Astronomy library
 *whereby the output from that function is transformed to be in terms of
 *hours minutes and seconds rather than as a fraction of a day.
 *The day from which the fraction is taken is assumed to have 24 hours,
 *with 60 minutes per hour and 60 seconds per minute.
 *
 *The UTC date is only pseudo-Julian, because there is not a fixed number
 *of seconds in a Julian day. The convention used in the IAU standard is
 *that the Julian day matches the UTC day regardless of whether the UTC day
 *is 86399, 86400 or 86401 SI seconds (depending on the presence of leap
 *seconds).
 *
 *Leap years are taken into account. February has 29 days if the year is
 *divisible by 4, unless the year is divisible by 100 unless the year is
 *divisible by 400. Leap seconds are also taken into account.
 *
 *UTC began at 1960 January 1.0 (JD 2436934.5) and this function should not
 *be called with an earlier date.
 *
 *The algorithm can be compiled for use in Matlab  using the
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[year,dayCount,second]=UTC2DayCount(Jul1,Jul2);
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"
#include "SharedTimeCode.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *Jul1,*Jul2;
    double *second;
    double *yearPtr, *dayCountPtr;
    size_t M, N, i,numElements;
    mxArray *yearMATLAB, *dayCountMATAB, *secondMATLAB;

    if(nrhs!=2){
        mexErrMsgTxt("Wrong number of inputs.");
    }
    
    if(nlhs>3) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
/*Extract the input arguments.*/
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    M=mxGetM(prhs[0]);
    N=mxGetN(prhs[0]);
    numElements=M*N;
    if(mxGetM(prhs[1])!=M||mxGetN(prhs[1])!=N) {
        mexErrMsgTxt("The dimensions of the inputs are inconsistent.");
    }
    
    Jul1 = (double*)mxGetData(prhs[0]);
    Jul2 = (double*)mxGetData(prhs[1]);
    
    //Allocate the return variables
    yearMATLAB=mxCreateDoubleMatrix(M,N,mxREAL);
    dayCountMATAB=mxCreateDoubleMatrix(M,N,mxREAL);
    secondMATLAB=mxCreateDoubleMatrix(M,N,mxREAL);
    
    yearPtr=(double*)mxGetData(yearMATLAB);
    dayCountPtr=(double*)mxGetData(dayCountMATAB);
    second=(double*)mxGetData(secondMATLAB);
    
    for(i=0;i<numElements;i++) {
        int retVal, year, dayCount;
        
        retVal=UTC2DayCountC(Jul1[i], Jul2[i], &year, &dayCount, &second[i]);
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
        yearPtr[i]=(double)year;
        dayCountPtr[i]=(double)dayCount;   
    }
    
    //Matlab users expect to be able to perform double precision arithmetic
    //on the results, so convert the results to doubles.
    plhs[0]=yearMATLAB;
    if(nlhs>1) {
        plhs[1]=dayCountMATAB;
        if(nlhs>2) {
            plhs[2]=secondMATLAB;
        } else {
            mxDestroyArray(secondMATLAB);
        }
    } else {
        mxDestroyArray(dayCountMATAB);
        mxDestroyArray(secondMATLAB);
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
