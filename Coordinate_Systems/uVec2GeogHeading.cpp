/**UVEC2GEOGHEADING  Obtain the geographic headings in radians East of true
*                 North as well as the elevations above the local tangent
*                 plane according to a particular reference ellipsoid that
*                 correspond to directions specified by unit vectors in
*                 ECEF coordinates.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*[geoEastOfNorth,angUpFromLevel]=uVec2GeogHeading(point,u);
*
*March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "CoordFuncs.hpp"
        
void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    
    if(mxGetNumberOfElements(prhs[0])<2||mxGetNumberOfElements(prhs[0])>3) {
        mexErrMsgTxt("point has the wrong dimensionality.");
        return;
    }
    
    const double *point=mxGetDoubles(prhs[0]);
    
    if(mxGetM(prhs[1])!=3) {
        mexErrMsgTxt("u has the wrong number of rows.");
        return;
    }
    
    const size_t numPoints=mxGetN(prhs[1]);
    const double *u=mxGetDoubles(prhs[1]);
    
    mxArray *geoEastOfNorthMATLAB=mxCreateDoubleMatrix(numPoints, 1, mxREAL);
    mxArray *angUpFromLevelMATLAB=mxCreateDoubleMatrix(numPoints, 1, mxREAL);
    double *geoEastOfNorth=mxGetDoubles(geoEastOfNorthMATLAB);
    double *angUpFromLevel=mxGetDoubles(angUpFromLevelMATLAB);

    uVec2GeogHeading(point, numPoints, u, geoEastOfNorth, angUpFromLevel);
    plhs[0]=geoEastOfNorthMATLAB;
    if(nlhs>1) {
        plhs[1]=angUpFromLevelMATLAB;
    } else {
        mxDestroyArray(angUpFromLevelMATLAB);
    }
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
