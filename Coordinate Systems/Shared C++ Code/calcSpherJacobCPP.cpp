/*CALCSPHERJACOBCPP  A C++ function to compute the Jacobian matrix for
 *                   spherical coordinates.
 *
 *INPUTS: J    A pointer to an array of doubles with 9 elements to hold the
 *             result. The result is stored by column. The first column
 *             has derivatives with respect to x, the second y and the
 *             third z. The rows correspond to spherical radius, azimuth
 *             and elevation. 
 *     point   The 3X1 point the format [range;azimuth;elevation], where
 *             the two angles are given in radians.
 *systemType   An integer specifying the axes from which the angles are
 *             measured. Possible vaues are
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
 *OUTPUTS: None. The results are placed in J.
 *
 *Further comments are given in the Matlab file calcSpherJacob.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
#include <math.h>

void calcSpherJacobCPP(double *J, const double *point, const int systemType) {
    double CartPoint[3],x,y,z,r,r2;
    
    spher2CartCPP(CartPoint,point,systemType);
    x=CartPoint[0];
    y=CartPoint[1];
    z=CartPoint[2];

    r=point[0];
    r2=r*r;
    
    if(systemType==0) {
        double x2y2Sum,xyDist;
        x2y2Sum=x*x+y*y;
        xyDist=sqrt(x2y2Sum);
        
        //Derivatives with respect to x.
        J[0]=x/r;
        J[1]=-y/x2y2Sum;
        J[2]=-x*z/(r2*xyDist);

        //Derivatives with respect to y.
        J[3]=y/r;
        J[4]=x/x2y2Sum;
        J[5]=-y*z/(r2*xyDist);

        //Derivatives with respect to z.
        J[6]=z/r;
        J[7]=0;
        J[8]=xyDist/r2;
    } else {
        double z2x2Sum,zxDist;
        z2x2Sum=z*z+x*x;
        zxDist=sqrt(z2x2Sum);
        
        //Derivatives with respect to x.
        J[0]=x/r;
        J[1]=z/z2x2Sum;
        J[2]=-x*y/(r2*zxDist);
        
        //Derivatives with respect to y.
        J[3]=y/r;
        J[4]=0;
        J[5]=zxDist/r2;
        
        //Derivatives with respect to z.
        J[6]=z/r;
        J[7]=-x/z2x2Sum;
        J[8]=-z*y/(r2*zxDist);
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
