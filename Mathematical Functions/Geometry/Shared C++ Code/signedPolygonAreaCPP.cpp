/*SIGNEDPOLYGONAREACPP A direct C++ implementation implementation of an
 *                  algorithm to calculate the signed area of a
 *                  non-self-intersecting 2D polygon given its vertices.
 *                  More comments on the implementation are given in the
 *                  Matlab implementation signedPolygonArea.m and its C++
 *                  version signedPolygonAreaCPP.cpp, which calls this
 *                  subroutine.
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathGeometricFuncs.hpp"

double signedPolygonAreaCPP(const double *vertices,const size_t numVertices) {
    double A;
    const double *prevVertex;
    size_t curIdx;
    
    A=0;
    prevVertex=vertices+2*(numVertices-1);
    for(curIdx=0;curIdx<numVertices;curIdx++) {
        const double *curVertex=vertices+2*curIdx;
        //A=A+det([prevVertex,curVertex]);
        A+=prevVertex[0]*curVertex[1]-prevVertex[1]*curVertex[0];
        prevVertex=curVertex;
    }
    A=A/2;
    return A;
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
