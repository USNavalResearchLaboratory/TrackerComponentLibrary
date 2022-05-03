/**GEOHEADING2UVEC Convert geographic headings in radians East of true
*             North with elevations above the local tangent plane at a
*             particular latitude and longitude with respect to a reference
*             ellipsoid to unit vectors in ECEF coordinates.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*u=geogHeading2uVec(point,geoEastOfNorth,angUpFromLevel)
*
*March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"
#include "CoordFuncs.hpp"
#include <algorithm>

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {    
    if(nrhs!=3) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    checkRealDoubleArray(prhs[2]);
    
    if(mxGetM(prhs[0])!=2&&mxGetM(prhs[0])!=3) {
        mexErrMsgTxt("point has the wrong dimensionality.");
        return;
    }
    const double *point=mxGetDoubles(prhs[0]);
    const size_t NPts=mxGetN(prhs[0]);
    const size_t NGeo=mxGetNumberOfElements(prhs[1]);
    const size_t NEl=mxGetNumberOfElements(prhs[2]);
    const double *geoEastOfNorth=mxGetDoubles(prhs[1]);
    const double *angUpFromLevel=mxGetDoubles(prhs[2]);
    
    if((NGeo!=NEl)&&(NGeo!=1)&&(NEl!=1)) {
        mexErrMsgTxt("The sizes of geoEastOfNorth and angUpFromLevel are inconsistent.");
        return;
    }

    if((NPts!=NEl)&&(NPts!=1)&&(NEl!=1)) {
        mexErrMsgTxt("The sizes of point and angUpFromLevel are inconsistent.");
        return;
    }

    if((NGeo!=NPts)&&(NGeo!=1)&&(NPts!=1)) {
        mexErrMsgTxt("The sizes of geoEastOfNorth and point are inconsistent.");
        return;
    }

    mxArray *uMATLAB=mxCreateDoubleMatrix(3,std::max(NPts,std::max(NGeo,NGeo)),mxREAL);
    double *u=mxGetDoubles(uMATLAB);

    geogHeading2uVecCPP(NPts,point, NGeo, geoEastOfNorth, NEl, angUpFromLevel, u);

    plhs[0]=uMATLAB;
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
