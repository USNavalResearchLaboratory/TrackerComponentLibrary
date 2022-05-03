/*A C++-only implementations of functions for computing the gradient of
*spherical azimuth and elevation angles. See the Matlab equivalent for
*more comments.
*
*February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
#include <math.h>

void spherAngGradientGenCPP(double *retMat,const double *xG,const size_t systemType,const double *lRx, const double *M) {
    double x,y,z;
    double temp[3];
    double J[6];
    //Transform from global coordinates to local coordinates.
    //xLocal=M*(xG(1:3)-lRx(1:3));
    temp[0]=xG[0]-lRx[0];
    temp[1]=xG[1]-lRx[1];
    temp[2]=xG[2]-lRx[2];
    x=M[0]*temp[0]+M[3]*temp[1]+M[6]*temp[2];
    y=M[1]*temp[0]+M[4]*temp[1]+M[7]*temp[2];
    z=M[2]*temp[0]+M[5]*temp[1]+M[8]*temp[2];

    if(systemType==0) {
        double r2=x*x+y*y+z*z;
        double sqrVal=x*x+y*y;
        double sqrtVal=sqrt(sqrVal);
        double denom=r2*sqrtVal;
        //Derivatives with respect to x.
        J[0]=-y/sqrVal;
        J[1]=-x*z/denom;

        //Derivatives with respect to y.
        J[2]=x/sqrVal;
        J[3]=-y*z/denom;

        //Derivatives with respect to z.
        J[4]=0;
        J[5]=sqrtVal/r2;
    } else if(systemType==2) {
        double r2=x*x+y*y+z*z;
        double sqrVal=x*x+y*y;
        double sqrtVal=sqrt(sqrVal);
        double denom=r2*sqrtVal;
        //Derivatives with respect to x.
        J[0]=-y/sqrVal;
        J[1]=x*z/denom;

        //Derivatives with respect to y.
        J[2]=x/sqrVal;
        J[3]=y*z/denom;

        //Derivatives with respect to z.
        J[4]=0;
        J[5]=-sqrtVal/r2;
    } else if (systemType==3) {
        double r2=x*x+y*y+z*z;
        double sqrVal=x*x+y*y;
        double sqrtVal=sqrt(sqrVal);
        double denom=r2*sqrtVal;
        //Derivatives with respect to x.
        J[0]=y/sqrVal;
        J[1]=-x*z/denom;

        //Derivatives with respect to y.
        J[2]=-x/sqrVal;
        J[3]=-y*z/denom;

        //Derivatives with respect to z.
        J[4]=0;
        J[5]=sqrtVal/r2;        
    }else{//Assume systemType==1
        double r2=x*x+y*y+z*z;
        double sqrVal=z*z+x*x;
        double sqrtVal=sqrt(sqrVal);
        double denom=r2*sqrtVal;
        //Derivatives with respect to x.
        J[0]=z/sqrVal;
        J[1]=-x*y/denom;

        //Derivatives with respect to y.
        J[2]=0;
        J[3]=sqrtVal/r2;

        //Derivatives with respect to z.
        J[4]=-x/sqrVal;
        J[5]=-z*y/denom;
    }
    
    //Rotate from local back to global coordinates.
    //J=J*M;
    retMat[0]=J[0]*M[0]+J[2]*M[1]+J[4]*M[2];
    retMat[1]=J[1]*M[0]+J[3]*M[1]+J[5]*M[2];
    retMat[2]=J[0]*M[3]+J[2]*M[4]+J[4]*M[5];
    retMat[3]=J[1]*M[3]+J[3]*M[4]+J[5]*M[5];
    retMat[4]=J[0]*M[6]+J[2]*M[7]+J[4]*M[8];
    retMat[5]=J[1]*M[6]+J[3]*M[7]+J[5]*M[8];
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
