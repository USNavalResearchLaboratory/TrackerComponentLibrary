/*SPHER2CART  A C++ function to convert a point from sperical coordinates
 *            to Cartesian coordinates.
 *
 *INPUTS: cartPoint    A pointer to an array of doubles with 3 elements to
 *                     hold the result in [x,y,z] order.
 *        point        The 3X1 point in spherical coordinates to be
 *                     converted, ordered [range;azimuth;elevation].
 *      systemType     An integer specifying the axes from which
 *                     the angles are measured. Possible vaues are
 *                   0 Azimuth is measured counterclockwise from the x-axis
 *                     in the x-y plane. Elevation is measured up from the
 *                     x-y plane (towards the z-axis). This is consistent
 *                     with common spherical coordinate systems for
 *                     specifying longitude (azimuth) and geocentric
 *                     latitude (elevation).
 *                   1 Azimuth is measured counterclockwise from the z-axis
 *                     in the z-x plane. Elevation is measured up from the
 *                     z-x plane (towards the y-axis). This is consistent
 *                     with some spherical coordinate systems that use the
 *                     z-axis as the boresight direction of the radar.
 *
 *OUTPUTS: None. The results are placed in cartPoint.
 *
 *Further comments are given in the Matlab file spher2Cart.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
//For sin and cos
#include <math.h>

void spher2CartCPP(double *cartPoint,const double *point, const int systemType) {
    double r, azimuth, elevation;
    double cosEl;
    r=point[0];
    azimuth=point[1];
    elevation=point[2];
    
    cosEl=cos(elevation);

    if(systemType==0) {
        cartPoint[0]=r*cos(azimuth)*cosEl;
        cartPoint[1]=r*sin(azimuth)*cosEl;
        cartPoint[2]=r*sin(elevation);
    } else {
        cartPoint[0]=r*sin(azimuth)*cosEl;
        cartPoint[1]=r*sin(elevation);
        cartPoint[2]=r*cos(azimuth)*cosEl;
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
