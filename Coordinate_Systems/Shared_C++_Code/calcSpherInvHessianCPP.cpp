/**CALCSPHERINVHESSIANCPP  A C++-only implementations of a function for
 *          computing the Hessian of of a 3D Cartesian point with respect
 *          to spherical azimuth and elevation. See the Matlab equivalent
 *          for more comments.
 *
 *July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//Needed for sin, cos
#include <cmath>
#include "CoordFuncs.hpp"

void calcSpherInvHessianCPP(double *H,const double *z,const size_t systemType) {
    const double r=z[0];
    const double sinAz=sin(z[1]);
    const double cosAz=cos(z[1]);
    const double sinEl=sin(z[2]);
    const double cosEl=cos(z[2]);

    if(systemType==0) {
        //d2xdrdr, H(1,1,1)
        H[0]=0;
        //d2xdAzdAz, H(2,2,1)
        H[4]=-r*cosAz*cosEl;
        //d2xdEldEl, H(3,3,1)
        H[8]=-r*cosAz*cosEl;
        //d2xdrdAz, H(1,2,1)
        H[3]=-cosEl*sinAz;
        //d2xdAzdr, H(2,1,1)
        H[1]=H[3];
        //d2xdrdEl, H(1,3,1)
        H[6]=-cosAz*sinEl;
        //d2xdEldr, H(3,1,1)
        H[2]=H[6];
        //d2xdAzdEl, H(2,3,1)
        H[7]=r*sinAz*sinEl;
        //d2xdEldAz, H(3,2,1)
        H[5]=H[7];

        //d2ydrdr, H(1,1,2)
        H[9]=0;
        //d2ydAzdAz, H(2,2,2)
        H[13]=-r*cosEl*sinAz;
        //d2ydEldEl, H(3,3,2)
        H[17]=-r*cosEl*sinAz;
        //d2ydrdAz, H(1,2,2)
        H[12]=cosAz*cosEl;
        //d2ydAzdr, H(2,1,2)
        H[10]=H[12];
        //d2ydrdEl, H(1,3,2)
        H[15]=-sinAz*sinEl;
        //d2ydEldr, H(3,1,2)
        H[11]=H[15];
        //d2ydAzdEl, H(2,3,2)
        H[16]=-r*cosAz*sinEl;
        //d2ydEldAz, H(3,2,2)
        H[14]=H[16];

        //d2zdrdr, H(1,1,3)
        H[18]=0;
        //d2zdAzdAz, H(2,2,3)
        H[22]=0;
        //d2zdEldEl, H(3,3,3)
        H[26]=-r*sinEl;
        //d2zdrdAz, H(1,2,3)
        H[21]=0;
        //d2zdAzdr, H(2,1,3)
        H[19]=H[21];
        //d2zdrdEl, H(1,3,3)
        H[24]=cosEl;
        //d2zdEldr, H(3,1,3)
        H[20]=H[24];
        //d2zdAzdEl, H(2,3,3)
        H[25]=0;
        //d2zdEldAz, H(3,2,3)
        H[23]=H[25];
    } else if(systemType==1) {
        //d2xdrdr, H(1,1,1)
        H[0]=0;
        //d2xdAzdAz, H(2,2,1)
        H[4]=-r*cosEl*sinAz;
        //d2xdEldEl, H(3,3,1)
        H[8]=-r*cosEl*sinAz;
        //d2xdrdAz, H(1,2,1)
        H[3]=cosAz*cosEl;
        //d2xdAzdr, H(2,1,1)
        H[1]=H[3];
        //d2xdrdEl, H(1,3,1)
        H[6]=-sinAz*sinEl;
        //d2xdEldr, H(3,1,1)
        H[2]=H[6];
        //d2xdAzdEl, H(2,3,1)
        H[7]=-r*cosAz*sinEl;
        //d2xdEldAz, H(3,2,1)
        H[5]=H[7];

        //d2ydrdr, H(1,1,2)
        H[9]=0;
        //d2ydAzdAz, H(2,2,2)
        H[13]=0;
        //d2ydEldEl, H(3,3,2)
        H[17]=-r*sinEl;
        //d2ydrdAz, H(1,2,2)
        H[12]=0;
        //d2ydAzdr, H(2,1,2)
        H[10]=H[12];
        //d2ydrdEl, H(1,3,2)
        H[15]=cosEl;
        //d2ydEldr, H(3,1,2)
        H[11]=H[15];
        //d2ydAzdEl, H(2,3,2)
        H[16]=0;
        //d2ydEldAz, H(3,2,2
        H[14]=H[16];

        //d2zdrdr, H(1,1,3)
        H[18]=0;
        //d2zdAzdAz, H(2,2,3)
        H[22]=-r*cosAz*cosEl;
        //d2zdEldEl, H(3,3,3)
        H[26]=-r*cosAz*cosEl;
        //d2zdrdAz, H(1,2,3)
        H[21]=-cosEl*sinAz;
        //d2zdAzdr, H(2,1,3)
        H[19]=H[21];
        //d2zdrdEl, H(1,3,3)
        H[24]=-cosAz*sinEl;
        //d2zdEldr, H(3,1,3)
        H[20]=H[24];
        //d2zdAzdEl, H(2,3,3)
        H[25]=r*sinAz*sinEl;
        //d2zdEldAz, H(3,2,3)
        H[23]=H[25];
    } else if(systemType==3) {
        //d2xdrdr, H(1,1,1)
        H[0]=0;
        //d2xdAzdAz, H(2,2,1)
        H[4]=-r*sinAz*cosEl;
        //d2xdEldEl, H(3,3,1)
        H[8]=-r*sinAz*cosEl;
        //d2xdrdAz, H(1,2,1)
        H[3]=cosEl*cosAz;
        //d2xdAzdr, H(2,1,1)
        H[1]=H[3];
        //d2xdrdEl, H(1,3,1)
        H[6]=-sinAz*sinEl;
        //d2xdEldr, H(3,1,1)
        H[2]=H[6];
        //d2xdAzdEl, H(2,3,1)
        H[7]=-r*cosAz*sinEl;
        //d2xdEldAz, H(3,2,1)
        H[5]=H[7];

        //d2ydrdr, H(1,1,2)
        H[9]=0;
        //d2ydAzdAz, H(2,2,2)
        H[13]=-r*cosEl*cosAz;
        //d2ydEldEl, H(3,3,2)
        H[17]=-r*cosEl*cosAz;
        //d2ydrdAz, H(1,2,2)
        H[12]=-sinAz*cosEl;
        //d2ydAzdr, H(2,1,2)
        H[10]=H[12];
        //d2ydrdEl, H(1,3,2)
        H[15]=-cosAz*sinEl;
        //d2ydEldr, H(3,1,2)
        H[11]=H[15];
        //d2ydAzdEl, H(2,3,2)
        H[16]=r*sinAz*sinEl;
        //d2ydEldAz, H(3,2,2)
        H[14]=H[16];

        //d2zdrdr, H(1,1,3)
        H[18]=0;
        //d2zdAzdAz, H(2,2,3)
        H[22]=0;
        //d2zdEldEl, H(3,3,3)
        H[26]=-r*sinEl;
        //d2zdrdAz, H(1,2,3)
        H[21]=0;
        //d2zdAzdr, H(2,1,3)
        H[19]=H[21];
        //d2zdrdEl, H(1,3,3)
        H[24]=cosEl;
        //d2zdEldr, H(3,1,3)
        H[20]=H[24];
        //d2zdAzdEl, H(2,3,3)
        H[25]=0;
        //d2zdEldAz, H(3,2,3)
        H[23]=H[25]; 
    } else {//systemType==2
        //d2xdrdr, H(1,1,1)
        H[0]=0;
        //d2xdAzdAz, H(2,2,1)
        H[4]=-r*cosAz*sinEl;
        //d2xdEldEl, H(3,3,1)
        H[8]=-r*cosAz*sinEl;
        //d2xdrdAz, H(1,2,1)
        H[3]=-sinAz*sinEl;
        //d2xdAzdr, H(2,1,1)
        H[1]=H[3];
        //d2xdrdEl, H(1,3,1)
        H[6]=cosAz*cosEl;
        //d2xdEldr, H(3,1,1)
        H[2]=H[6];
        //d2xdAzdEl, H(2,3,1)
        H[7]=-r*cosEl*sinAz;
        //d2xdEldAz, H(3,2,1)
        H[5]=H[7];

        //d2ydrdr, H(1,1,2)
        H[9]=0;
        //d2ydAzdAz, H(2,2,2)
        H[13]=-r*sinAz*sinEl;
        //d2ydEldEl, H(3,3,2)
        H[17]=-r*sinAz*sinEl;
        //d2ydrdAz, H(1,2,2)
        H[12]=cosAz*sinEl;
        //d2ydAzdr, H(2,1,2)
        H[10]=H[12];
        //d2ydrdEl, H(1,3,2)
        H[15]=cosEl*sinAz;
        //d2ydEldr, H(3,1,2)
        H[11]=H[15];
        //d2ydAzdEl, H(2,3,2)
        H[16]=r*cosAz*cosEl;
        //d2ydEldAz, H(3,2,2)
        H[14]=H[16];

        //d2zdrdr, H(1,1,3)
        H[18]=0;
        //d2zdAzdAz, H(2,2,3)
        H[22]=0;
        //d2zdEldEl, H(3,3,3)
        H[26]=-r*cosEl;
        //d2zdrdAz, H(1,2,3)
        H[21]=0;
        //d2zdAzdr, H(2,1,3)
        H[19]=H[21];
        //d2zdrdEl, H(1,3,3)
        H[24]=-sinEl;
        //d2zdEldr, H(3,1,3)
        H[20]=H[24];
        //d2zdAzdEl, H(2,3,3
        H[25]=0;
        //d2zdEldAz, H(3,2,3)
        H[23]=H[25];
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
