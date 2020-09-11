/*TWOLINEINTERSECTIONPOINT2D A C++ implementation of a function that given
 *                    two (infinite) 2D lines that are not parallel, finds
 *                    the point of intersection. A line is specified by
 *                    providing two points on the line.
 *
 *INPUTS: line1 A 2X2 matrix of two points in the first line where
 *              line1(1,:) are the x-coordinates and line1(2,:) are the
 *               y-coordinates.
 *         line2 A 2X2 matrix of two points in the second line defined the
 *               same way as line1.
 *
 *OUTPUTS: point The 2X1 intersection point of the two lines given as [x;y]
 *               components.
 *
 *The formula in terms of matrix determinants is taken from [1].
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *point=twoLineIntersectionPoint2D(line1,line2);
 *
 *REFERENCES:
 *[1] Weisstein, Eric W. "Line-Line Intersection." From MathWorld--A
 *    Wolfram Web Resource.
 *    http://mathworld.wolfram.com/Line-LineIntersection.html
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//This header is required by Matlab.
#include "mex.h"
//This is for input validation
#include "MexValidation.h"
#include "mathGeometricFuncs.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *line1,*line2,*point;
    mxArray *pointMat;
    
    if(nrhs!=2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Incorrect number of outputs.");
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);

    if(mxGetM(prhs[0])!=2||mxGetN(prhs[0])!=2||mxGetM(prhs[1])!=2||mxGetN(prhs[1])!=2) {
        mexErrMsgTxt("The inputs are not the correct size.");
        return;  
    }
    
    line1=mxGetDoubles(prhs[0]);
    line2=mxGetDoubles(prhs[1]);

    //Allocate space for the return point.
    pointMat=mxCreateDoubleMatrix(2,1,mxREAL);
    point=mxGetDoubles(pointMat);
    
    twoLineIntersectionPoint2DCPP(line1,line2,point);
    
    //Set the return value.
    plhs[0]=pointMat;
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
