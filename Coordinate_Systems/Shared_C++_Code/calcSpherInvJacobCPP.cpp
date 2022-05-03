/**CALCSPHERINVJACOBCPP  A C++-only implementations of a function for
 *          computing the Jacobian of of a 3D Cartesian point with respect
 *          to spherical azimuth and elevation. See the Matlab equivalent
 *          for more comments.
 *
 *July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//Needed for sin, cos
#include <cmath>
#include "CoordFuncs.hpp"

void calcSpherInvJacobCPP(double *J,const double *z,const size_t systemType) {
    const double r=z[0];
    const double sinAz=sin(z[1]);
    const double cosAz=cos(z[1]);
    const double sinEl=sin(z[2]);
    const double cosEl=cos(z[2]);
    
    if(systemType==0) {
        //dx/dr, J(1,1)
        J[0]=cosAz*cosEl;
        //dy/dr, J(2,1)
        J[1]=cosEl*sinAz;
        //dz/dr, J(3,1)
        J[2]=sinEl;
        //dx/dAz, J(1,2)
        J[3]=-r*cosEl*sinAz;
        //dy/dAz, J(2,2)
        J[4]=r*cosAz*cosEl;
        //dz/dAz, J(3,2)
        J[5]=0;
        //dx/dEl, J(1,3)
        J[6]=-r*cosAz*sinEl;
        //dy/dEl, J(2,3)
        J[7]=-r*sinAz*sinEl;
        //dz/dEl, J(3,3)
        J[8]=r*cosEl;
    } else if(systemType==1) {
        //dx/dr, J(1,1)
        J[0]=cosEl*sinAz;
        //dy/dr, J(2,1)
        J[1]=sinEl;
        //dz/dr, J(3,1)
        J[2]=cosAz*cosEl;
        //dx/dAz, J(1,2)
        J[3]=r*cosAz*cosEl;
        //dy/dAz, J(2,2)
        J[4]=0;
        //dz/dAz, J(3,2)
        J[5]=-r*cosEl*sinAz;
        //dx/dEl, J(1,3)
        J[6]=-r*sinAz*sinEl;
        //dy/dEl, J(2,3)
        J[7]=r*cosEl;
        //dz/dEl, J(3,3)
        J[8]=-r*cosAz*sinEl;
    } else if (systemType==3) {
        //dx/dr, J(1,1)
        J[0]=sinAz*cosEl;
        //dy/dr, J(2,1)
        J[1]=cosAz*cosEl;
        //dz/dr, J(3,1)
        J[2]=sinEl;
        //dx/dAz, J(1,2)
        J[3]=r*cosEl*cosAz;
        //dy/dAz, J(2,2)
        J[4]=-r*sinAz*cosEl;
        //dz/dAz, J(3,2)
        J[5]=0;
        //dx/dEl, J(1,3)
        J[6]=-r*sinAz*sinEl;
        //dy/dEl, J(2,3)
        J[7]=-r*cosAz*sinEl;
        //dz/dEl, J(3,3)
        J[8]=r*cosEl;
    } else {//systemType==2
        //dx/dr, J(1,1)
        J[0]=cosAz*sinEl;
        //dy/dr, J(2,1)
        J[1]=sinAz*sinEl;
        //dz/dr, J(3,1)
        J[2]=cosEl;
        //dx/dAz, J(1,2)
        J[3]=-r*sinAz*sinEl;
        //dy/dAz, J(2,2)
        J[4]=r*cosAz*sinEl;
        //dz/dAz, J(3,2)
        J[5]=0;
        //dx/dEl, J(1,3)
        J[6]=r*cosAz*cosEl;
        //dy/dEl, J(2,3)
        J[7]=r*cosEl*sinAz;
        //dz/dEl, J(3,3)
        J[8]=-r*sinEl;
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
