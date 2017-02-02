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
    //detL1=det(line1);
    const double detL1=line1[0]*line1[3]-line1[1]*line1[2];
    //detL2=det(line2);
    const double detL2=line2[0]*line2[3]-line2[1]*line2[2];
    //deltaL1=-diff(line1,1,2);
    const double deltaL1[2]={line1[0]-line1[2],line1[1]-line1[3]};
    //deltaL2=-diff(line2,1,2);
    const double deltaL2[2]={line2[0]-line2[2],line2[1]-line2[3]};
    //denom=det([deltaL1,deltaL2]);
    const double denom=deltaL1[0]*deltaL2[1]-deltaL1[1]*deltaL2[0];
    
    //x=det([[detL1;detL2],[deltaL1(1);deltaL2(1)]])/denom;
    point[0]=(detL1*deltaL2[0]-detL2*deltaL1[0])/denom;
    //y=det([[detL1;detL2],[deltaL1(2);deltaL2(2)]])/denom;
    point[1]=(detL1*deltaL2[1]-detL2*deltaL1[1])/denom;
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
