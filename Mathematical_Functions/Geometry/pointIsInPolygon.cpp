/*POINTISINPOLYGON A C++ implementation of a function that given a polygon
*               specified by a number of vertices, determines whether a
*               point is in the polygon. A simple polygon has lines that
*               do not cross. If the polygon is not simple (is
*               self-intersecting), then the even-odd rule is used to
*               determine whether a point is in the polygon or not. As the
*               algorithm is a type of winding number algorithm, the
*               winding number omega can also be obtained.
*
*INPUTS: vertices A 2XN matrix of the N vertices making up the polygon in
*                 order. Edges are between neighboring vertices. The last
*                 vertex can be the same as the first vertex. If not, it is
*                 assumed that an edge exists between the last and first
*                 vertices.
*           point A 2XnumPoints set of points that will be determined to be
*                 inside or outside of the polygon given by vertices.
* boundaryIsImportant An optional boolean variable indicating whether the
*                 boundary of the polygon is important. If
*                 boundaryIsImportant=true, then an algorithm that will
*                 correctly indicate points on the boundary as being in
*                 the polygon will be used. If it is false, then the
*                 results for points on the boundary can be inconsistent,
*                 though the algorithm will be slightly faster. The default
*                 is true.
*
*OUTPUTS: isInPolygon An NX1 vector where the ith element is true if the
*                     ith is true if the ith point is in the polygon, false
*                     otherwise. If boundaryIsImportant=false, the results
*                     can be inconsistent on the boundary of the polygon.
*              omegas An NX1 vector of the integer winding number obtained
*                     for each point. This is not meainingful for points on
*                     the boundary of the polygon. 
*
*The implementation is Algorithms 6 and 7 from [1].
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
* The algorithm is run in Matlab using the command format
* [isInPolygon,omegas]=pointIsInPolygon(vertices,point,boundaryIsImportant);
* or
* [isInPolygon,omegas]=pointIsInPolygon(vertices,point);
*
*REFERENCES:
*[1] K. Hormann and A. Agathos, "The point in polygon problem for arbitrary
*    polygons," Computational Geometry, vol. 20, no. 3, pp. 131-144, Nov.
*    2001.
*
*December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "mathGeometricFuncs.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    bool boundaryIsImportant=true;
    size_t numRow, numCol, numVertices;
    size_t numPoints;
    mxArray *isInPolygonMatlab, *omegasMatlab;
    ptrdiff_t *omegas;
    mxLogical *isInPolygon;
    double *P, *R;
    size_t i;
    
    if(nrhs<2){
        mexErrMsgTxt("Not enough inputs.");
        return;
    }
    
    if(nrhs>3) {
        mexErrMsgTxt("Too many inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Too many outputs.");
    }

    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    P=mxGetDoubles(prhs[0]);//The vertex data
    R=mxGetDoubles(prhs[1]);//The points
    
    numRow=mxGetM(prhs[0]);
    numCol=mxGetN(prhs[0]);

    if(numRow !=2) {
        mexErrMsgTxt("The vertices must be two-dimensional.");
    }

    if(numCol<3) {
        mexErrMsgTxt("There must be at least 3 vertices.");
    }
    numVertices=numCol;
    
    numRow=mxGetM(prhs[1]);
    numPoints=mxGetN(prhs[1]);
    if(numRow!=2||numPoints==0) {
        mexErrMsgTxt("The points must be two-dimensional.");
    }
    
    if(nrhs>2) {
        boundaryIsImportant=getBoolFromMatlab(prhs[2]);
    }
    
    //Allocate space for the return values.
    omegasMatlab=allocSignedSizeMatInMatlab(numPoints,1);
    isInPolygonMatlab=mxCreateLogicalMatrix(numPoints,1);
    isInPolygon=mxGetLogicals(isInPolygonMatlab);
    
    if(sizeof(ptrdiff_t)==4) {//32 bit
        omegas=reinterpret_cast<ptrdiff_t*>(mxGetInt32s(omegasMatlab));
    } else {//64 bit
        omegas=reinterpret_cast<ptrdiff_t*>(mxGetInt64s(omegasMatlab));
    }

    for(i=0;i<numPoints;i++) {
        isInPolygon[i]=static_cast<mxLogical>(pointIsInPolygonCPP(P,numVertices,R+2*i,boundaryIsImportant,omegas+i));
    }

    //Set the return values
    plhs[0]=isInPolygonMatlab;
    
    if(nlhs>1) {
        plhs[1]=omegasMatlab;
    } else {
        mxDestroyArray(omegasMatlab);
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
