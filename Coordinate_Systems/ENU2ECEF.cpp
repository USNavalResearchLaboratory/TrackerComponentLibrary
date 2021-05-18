/**ENU2ECEF Convert from a local East-North-Up (ENU) Cartesian cooridnate
*          system to an Earth-centered Earth-fixed (ECEF) Cartesian
*          coordinate system. See the Matlab implementation for more
*          details.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*xECEF=ENU2ECEF(plhOrigin,xENU,a,f)
*
*March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"
#include "CoordFuncs.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double a, f;
    mxArray *uMATLAB,*MMATLAB=NULL;
    double *u, *M;
    double MLoc[9];

    if(nrhs<2||nrhs>4) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        a=getDoubleFromMatlab(prhs[2]);
    } else {//Load the default value if none is supplied.
        a=getScalarMatlabClassConst("Constants","WGS84SemiMajorAxis");
    }

    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        f=getDoubleFromMatlab(prhs[3]);
    } else {//Load the default value if none is supplied.
        f=getScalarMatlabClassConst("Constants","WGS84Flattening");
    }

    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    
    if(mxGetNumberOfElements(prhs[0])!=3) {
        mexErrMsgTxt("plhOrigin has the wrong dimensionality.");
        return;
    }
    
    const double *plhOrigin=mxGetDoubles(prhs[0]);
    
    if(mxGetM(prhs[1])!=3) {
        mexErrMsgTxt("xENU has the wrong number of rows.");
        return;
    }
    
    const size_t numPts=mxGetN(prhs[1]);
    const double *xENU=mxGetDoubles(prhs[1]);
    
    uMATLAB=mxCreateDoubleMatrix(3, numPts,mxREAL);
    u=mxGetDoubles(uMATLAB);
    plhs[0]=uMATLAB;
    
    if(nlhs>1) {
        MMATLAB=mxCreateDoubleMatrix(3, 3,mxREAL);
        M=mxGetDoubles(MMATLAB);
        plhs[1]=MMATLAB;
    } else {
        M=MLoc;
    }
    
    ENU2ECEFCPP(numPts, xENU, plhOrigin, a, f, u, M);
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
