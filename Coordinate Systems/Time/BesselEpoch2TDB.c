/**BESSELEPOCH2TDB Convert a Besselian epoch into a two-part Julian date in
*                  barycentric dynamical time (TDB). Besselian years are
*                  related to the apparent position of the Sun and are
*                  generaly considered outdated.
*
*INPUTS: bessEpoch A matrix of Besselian epochs that are to be converted to
*                  Julian dates in TDB.
*
*OUTPUTS: Jul1, Jul2 Matrices of two parts of a Julian date given
*                    in TDB. The units of the date are days. The full
*                    date is the sum of both terms. The entries in the
*                    matrices match the corresponding entries in bessEpoch.
*
*The relationship between Julian dates in TDB and Besselian epochs is given
*in
*[1] J. H. Lieske, "Precession matrix based on IAU (1976) system of
*astronomical constants," Astronomy and Astrophysics, vol. 73, no. 3, pp.
*282-284, Mar. 1979.
*A Besselian epoch is a factional year number denominated in terms of a
*tropical year. For example, 1900.5 is a time near the middle of the year
*1900. A Besselian year is not the same duration as a Julian year in TT.
*
*This is a mex wrapper for the function iauEpb2jd in the International
*Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
*Though the IAU's library does not explicitely say that the Julian date
*must be in TDB, the Besselian epoch does not make sense unless the Julian
*date is in TDB as one can ascertain from [1]. 
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*[Jul1, Jul2]=BesselEpoch2TDB(bessEpoch);
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
    double *BesselEpoch,*Jul1Ret,*Jul2Ret;
    size_t numRow, numCol, numElements,i;
    mxArray *Jul1RetMATLAB,*Jul2RetMATLAB;
    
    if(nrhs!=1){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>2){
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    numRow=mxGetM(prhs[0]);
    numCol=mxGetN(prhs[0]);
    numElements=numRow*numCol;
    
    if(numElements==0) {
        mexErrMsgTxt("The dimensionalities of the inputs are incorrect.");
        return;
    }

    checkRealDoubleArray(prhs[0]);
    BesselEpoch=(double*)mxGetData(prhs[0]);

    //Allocate space for the return variables.
    Jul1RetMATLAB=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    Jul1Ret=(double*)mxGetData(Jul1RetMATLAB);
    Jul2RetMATLAB=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    Jul2Ret=(double*)mxGetData(Jul2RetMATLAB);

    for(i=0;i<numElements;i++) {
        /*Call the function in the SOFA library.*/
        iauEpb2jd(BesselEpoch[i], Jul1Ret+i, Jul2Ret+i);
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
