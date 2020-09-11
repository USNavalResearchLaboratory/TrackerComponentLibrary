/*TWOLINEINTERSECTIONPOINT2DCPP A direct C++ implementation of a function
 *                that given two (infinite) 2D lines that are not parallel,
 *                finds the point of intersection. A line is specified by
 *                providing two points on the line. More details on the
 *                algorithms are given in the C++ interface for Matlab
 *                twoLineIntersectionPoint2D.cpp and the Matlab function
 *                twoLineIntersectionPoint2D.m
 *
 *December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "mathGeometricFuncs.hpp"

void twoLineIntersectionPoint2DCPP(const double *line1, const double *line2,double *point) {
    const double x1=line1[0];
    const double y1=line1[1];
    const double x2=line1[2];
    const double y2=line1[3];
    const double x3=line2[0];
    const double y3=line2[1];
    const double x4=line2[2];
    const double y4=line2[3];
    double num, den;
    
    num=(x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4); 
    den=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
    point[0]=num/den;//x

    num=(x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4); 
    den=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4); 
    point[1]=num/den;//y
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
