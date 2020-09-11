/*SIGNEDPOLYGONAREA A C++ implementation of an algorithm to calculate the
 *                  signed area of a non-self-intersecting 2D polygon given
 *                  its vertices. The magnitude of the area is the typical
 *                  definition of area that one might expect. The sign of
 *                  the area is positive if the points are arranged in a
 *                  counterclockwise order around the polygon and negative
 *                  if the points are in a clockwise order.
 *
 *INPUTS: vertices A 2XN list of the N vertices of the polygon in order
 *                 around the polygon. Duplicate vertices are allowed. It
 *                 does not matter whether or not the first vertex is
 *                 repeated at the end.
 *
 *OUTPUTS:       A The signed area of the polygon.
 *
 *The formula for computing the signed area of a non-self-intersecting is
 *taken from [1[.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *A=signedPolygonArea(vertices);
 *
 *REFERENCES:
 *[1] Weisstein, Eric W. "Polygon Area." From MathWorld--A Wolfram Web 
 *    Resource. http://mathworld.wolfram.com/PolygonArea.html
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//This header is required by Matlab.
#include "mex.h"
//This is for input validation
#include "MexValidation.h"
//To determine the intersection point of two lines.
#include "mathGeometricFuncs.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t numVertices;
    double *vertices, A;

    if(nrhs!=1) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;
    }

    checkRealDoubleArray(prhs[0]);
    
    //If an empty matrix is passed, return zero area.
    if(mxIsEmpty(prhs[0])) {
        const double retVal=0;
        plhs[0]=doubleMat2Matlab(&retVal,1,1);
        return;
    }
    
    if(mxGetM(prhs[0])!=2) {
        mexErrMsgTxt("The points have the wrong dimensionality.");
        return;
    }
    
    numVertices=mxGetN(prhs[0]);
    vertices=mxGetDoubles(prhs[0]);
    
    A=signedPolygonAreaCPP(vertices,numVertices);
    plhs[0]=doubleMat2Matlab(&A,1,1);
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
