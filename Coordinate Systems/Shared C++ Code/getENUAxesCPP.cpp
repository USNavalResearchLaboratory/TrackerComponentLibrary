/*GETENUAXESCPP  A C++ function to compute the basis vectors of an
 *               East-North-Up coordinate system.
 *
 *INPUTS: u, c   Pointers to arrays in which the unit vectors for the ENU
 *               axes and the magnitudes of the unnormalized vectors are
 *               placed.
 *     plhPoint  The 3X1 array at which the axes are to be found given in
 *               terms of [latitude;longitude;height] with the geodetic 
 *               latitude and longitude in radians and the height in
 *               meters. The latitude should be between -pi/2 and pi/2.
 *  justVertical A boolean parameter. If true then u and c only for the Up
 *               direction will be returned.
 *           a   The semi-major axis of the reference ellipsoid.
 *           f   The flattening factor of the reference ellipsoid
 *
 *OUTPUTS: None. The results are placed in u and c.
 *
 *Further comments are given in the Matlab file getENUAxes.
 *
 *February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
#include <math.h>
#include <stddef.h>

void getENUAxesCPP(double *u, double *c,const double *plhPoint,const bool justVertical,const double a,const double f) {
    const double phi=plhPoint[0];//The latitude
    const double lambda=plhPoint[1];//The longitude
    const double h=plhPoint[2];//The height
    const double sinP=sin(phi);
    const double cosP=cos(phi);
    const double sinL=sin(lambda);
    const double cosL=cos(lambda);
    double *u1,*u2,*u3;//Pointer to the normalized East, North, and Up vectors
    double *c1,*c2,*c3;//Pointers to the magnitudes of the unnormalized ENU vectors.
    size_t i;
    
    if(justVertical==false) {
        u1=u;
        u2=u+3;
        u3=u+6;
        c1=c;
        c2=c+1;
        c3=c+2;
    } else {
        u3=u;
        c3=c;
        //These variables are not used when justVertical==true, but are
        //initialized to get rid of a warning when compiled with 
        //-Wconditional-uninitialized
        u1=NULL;
        u2=NULL;
        c1=NULL;
        c2=NULL;
    }
    
    //u3 is dr/dh (Up)
    u3[0]=cosP*cosL;
    u3[1]=cosP*sinL;
    u3[2]=sinP;
    *c3=sqrt(u3[0]*u3[0]+u3[1]*u3[1]+u3[2]*u3[2]);//Barring precision problems, this is always one.
    for(i=0;i<3;i++){u3[i]=u3[i]/(*c3);}
    
    if(justVertical==false) {
        //The square of the first numerical eccentricity
        const double e2=2*f-f*f;
        //The normal radius of curvature.
        const double Ne=a/sqrt(1-e2*sinP*sinP);
        //The derivative of the normal radius of curvature with respect to phi.
        const double dNedPhi=a*e2*cosP*sinP/pow(1-e2*sinP*sinP,3.0/2.0);

        //u1 is dr/dlambda, normalized (East).
        u1[0]=-(Ne+h)*cosP*sinL;
        u1[1]=(Ne+h)*cosP*cosL;
        u1[2]=0;
        *c1=sqrt(u1[0]*u1[0]+u1[1]*u1[1]+u1[2]*u1[2]);

        //u2 is dr/dphi, normalized (North)
        u2[0]=(cosP*dNedPhi-(Ne+h)*sinP)*cosL;
        u2[1]=(cosP*dNedPhi-(Ne+h)*sinP)*sinL;
        u2[2]=(Ne*(1-e2)+h)*cosP+(1-e2)*dNedPhi*sinP;
        *c2=sqrt(u2[0]*u2[0]+u2[1]*u2[1]+u2[2]*u2[2]);
        for(i=0;i<3;i++){u2[i]=u2[i]/(*c2);}

        //If the point is too close to the poles, then it is possible that c1 is
        //nearly equal to zero. However, u1 can just be found by orthogonality:
        //it is orthogonal to u3 and u2. u1=cross(u2,u3);
        u1[0]=u2[1]*u3[2]-u2[2]*u3[1];
        u1[1]=u2[2]*u3[0]-u2[0]*u3[2];
        u1[2]=u2[0]*u3[1]-u2[1]*u3[0];
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
