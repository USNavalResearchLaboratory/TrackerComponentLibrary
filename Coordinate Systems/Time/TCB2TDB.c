/**TCB2TDB Convert from  barycentric coordinate time (TCB) to barycentric
 *         dynamical time (TDB), both represented as two-part Julian dates.
 *
 *INPUTS: Jul1, Jul2 Two parts of a Julian date given in TCB. The units of
 *                   the date are days. The full date is the sum of both
 *                   terms. The date is broken into two parts to provide
 *                   more bits of precision. It does not matter how the
 *                   date is split.
 *
 *OUTPUTS: Jul1, Jul2 The time as a Julian date in TDB.
 *
 *This is a mex wrapper for the function iauTcbtdb in the International
 *Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[Jul1,Jul2]=TCB2TDB(Jul1,Jul2);
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
    double Jul1,Jul2;
    int retVal;

    if(nrhs!=2){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
/*Extract the input arguments.*/
    Jul1 = getDoubleFromMatlab(prhs[0]);
    Jul2 = getDoubleFromMatlab(prhs[1]);
    
/*Call the function in the SOFA library.*/
    retVal=iauTcbtdb(Jul1, Jul2, &Jul1, &Jul2);
    
    if(retVal!=0) {
        mexErrMsgTxt("An error occurred.");
    }
    
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *(double*)mxGetData(plhs[0])=Jul1;
    if(nlhs>1) {
        plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
        *(double*)mxGetData(plhs[1])=Jul2;
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
