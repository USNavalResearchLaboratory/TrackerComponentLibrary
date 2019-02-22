/**ELLIPS2CART A C-code function to convert ellipsoidal coordinates to ECEF
*              Cartesian coordinates.
*
*INPUTS: points One or more points given in geodetic latitude and
*               longitude, in radians, and height, in meters that are to be
*               converted to Cartesian coordinates. To convert N points,
*               points is a 3XN matrix with each column having the format
*               [latitude;longitude; height].
*             a The semi-major axis of the reference ellipsoid. If this
*               argument is omitted, the value in
*               Constants.WGS84SemiMajorAxis is used.
*             f The flattening factor of the reference ellipsoid. If this
*               argument is omitted, the value in Constants.WGS84Flattening
*               is used.
*
*OUTPUTS: cartPoints For N points, cartPoints is a 3XN matrix of the
*               converted points with each column having the format
*               [x;y;z].
*
*This just calls the appropriate conversion function from the International
*Astronomical Union's Standard's of Fundamental Astronomy library.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*cartPoints=ellips2Cart(points);
*or
*cartPoints=ellips2Cart(points,a);
*or
*cartPoints=ellips2Cart(points,a,f);
*
*Note that the WGS-84 ellipsoid parameters obtained from the IAU's library
*should match those in Constants.WGS84SemiMajorAxis and
*Constants.WGS84Flattening.
*
*December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t numRow,numVec,i;
    mxArray *retMat;
    double *points, a, f;
    double *retData;
    
    if(nrhs>3||nrhs<1){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    checkRealDoubleArray(prhs[0]);
    
    numRow = mxGetM(prhs[0]);
    numVec = mxGetN(prhs[0]);
    
    if(numRow!=3) {
        mexErrMsgTxt("The input vector has a bad dimensionality.");
    }
    
    points=(double*)mxGetData(prhs[0]);
    //points[0] is latitude
    //points[1] is longitude
    //points[2] is height.
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        a=getDoubleFromMatlab(prhs[1]);
    } else {
        a=getScalarMatlabClassConst("Constants", "WGS84SemiMajorAxis");
    }

    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        f=getDoubleFromMatlab(prhs[2]);
    } else {
        f=getScalarMatlabClassConst("Constants", "WGS84Flattening");   
    }
    
    retMat=mxCreateDoubleMatrix(3,numVec,mxREAL);
    retData=(double*)mxGetData(retMat);
    for(i=0;i<numVec;i++) {
        int retVal;
        retVal= iauGd2gce(a,f, points[3*i+1], points[3*i], points[2+3*i], retData+3*i);
        if(retVal!=0) {
            mexErrMsgTxt("A conversion error has occurred.");
        }
    }
    plhs[0]=retMat;
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

