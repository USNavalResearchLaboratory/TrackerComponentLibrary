/**MATHGEOMETRICFUNCS  A header file for C++ implementations of 
 *           mathematical functions related to computational geometry. 
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef MATHGEOMETRICFUNCSCPP
#define MATHGEOMETRICFUNCSCPP
#include <stddef.h>

/**POINTISINPOLYGONCPP Given a polygon specified by a number of vertices,
*              determine whether a point is in the polygon. A simple
*              polygon has lines that do not cross. If the polygon is not
*              simple (is self-intersecting), then the even-odd rule is
*              used to determine whether a point is in the polygon or
*              not. As the algorithm is a type of winding number algorithm,
*              the winding number omega is also be obtained.
*
*INPUTS: P An array holding 2*numVertices values, where each pair is a 2D
*          point in the polygon in order. The last vertex can be the same
*          as the first vertex. If not, it is assumed that an edge exists
*          between the last and first vertices.
* numVertices The number of vertices in P.
*        R A length 2 array holding the point that is to be tested for
*          being in the polygon.
* boundaryIsImportant An boolean value indicating whether the boundary of
*          the polygon is important. If boundaryIsImportant=true, then an
*          algorithm that will correctly indicate points on the boundary as
*          being in the polygon will be used. If it is false, then the
*          results for points on the boundary can be inconsistent, though
*          the algorithm will be slightly faster.
*    omega The winding number number obtained for the point. This is not
*          meaningful for points on the boundary of the polygon.
*
*OUTPUTS: The winding number is placed in omega. the return value is a
*         boolean value indicating whether the point is in the polygon.
*
*December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
bool pointIsInPolygonCPP(const double *P, const size_t numVertices, const double *R, const bool boundaryIsImportant,ptrdiff_t *omega);

/**TWOLINESINTERSECTIONPOINT2DCPP Given two (infinite) 2D lines that are
*                   not parallel, find the point of intersection. A line is
*                   specified by providing two points on the line.
*
*INPUTS: line1 A length-4 vector holding two (x,y) points, one after the
*              other, that are on the first line.
*        line2 A length-4 vector holding two (x,y) points, one after the
*              other, that are on the second line.
*        point A length-2 buffer into which the intersection point is
*              placed.
*
*OUTPUTS: The result is placed in point.
*
*December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
void twoLineIntersectionPoint2DCPP(const double *line1, const double *line2,double *point);

/**SIGNEDPOLYGONAREA Calculate the signed area of a non-self-intersecting
*                   2D polygon given its vertices. The magnitude of the
*                   area is the typical definition of area that one might
*                   expect. The sign of the area is positive if the points
*                   are arranged in a counterclockwise order around the
*                   polygon and negative if the points are in a clockwise
*                   order.
*
*INPUTS: vertices A length 2*numVertices buffer holding numVertices (x,y)
*                 vertices of the polygon in a particular order around the
*                 polygon. Duplicate vertices are allowed. It does not
*                 matter whether or not the first vertex is repeated at the
*                 end.
*    numVertices The number of vertices in vertices.
*
*OUTPUTS: The return value is the signed polygon area.
*
*The formula for computing the signed area of a non-self-intersecting
*polygon is taken from [1].
*
*REFERENCES:
*[1] Weisstein, Eric W. "Polygon Area." From MathWorld--A Wolfram Web
*    Resource. http://mathworld.wolfram.com/PolygonArea.html
*
*December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
double signedPolygonAreaCPP(const double *vertices,const size_t numVertices);

#endif

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
