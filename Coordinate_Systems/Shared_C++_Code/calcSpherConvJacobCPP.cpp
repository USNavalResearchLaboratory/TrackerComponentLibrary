/**CALCSPHERCONVJACOBCPP These are C++-only implementations of functions
 *for computing gradients with respect to certain parameters. See the
 *Matlab equivalents for more comments.
 *
*February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
#include <math.h>

void calcSpherConvJacobGenCPP(double *J,const double *zSpher,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M) {
double JRange[3];
double tempSpace[3];
double JAngle[6];
double zCart[3];

spher2CartGenCPP(zCart,zSpher,systemType,useHalfRange,lTx,lRx,M);
rangeGradientCPP(3,JRange,tempSpace,zCart,useHalfRange,lTx,lRx);
spherAngGradientGenCPP(JAngle,zCart,systemType,lRx,M);

J[0]=JRange[0];
J[1]=JAngle[0];
J[2]=JAngle[1];

J[3]=JRange[1];
J[4]=JAngle[2];
J[5]=JAngle[3];

J[6]=JRange[2];
J[7]=JAngle[4];
J[8]=JAngle[5];
}

void calcSpherConvJacobCPP(double *J, const double *point, const size_t systemType) {
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
    } else if(systemType==2) { 
        double x2y2Sum,xyDist;
        x2y2Sum=x*x+y*y;
        xyDist=sqrt(x2y2Sum);
        
        //Derivatives with respect to x.
        J[0]=x/r;
        J[1]=-y/x2y2Sum;
        J[2]=x*z/(r2*xyDist);

        //Derivatives with respect to y.
        J[3]=y/r;
        J[4]=x/x2y2Sum;
        J[5]=y*z/(r2*xyDist);

        //Derivatives with respect to z.
        J[6]=z/r;
        J[7]=0;
        J[8]=-xyDist/r2;
    } else if(systemType==3) {
        double x2y2Sum,xyDist;
        x2y2Sum=x*x+y*y;
        xyDist=sqrt(x2y2Sum);
        
        //Derivatives with respect to x.
        J[0]=x/r;
        J[1]=y/x2y2Sum;
        J[2]=-x*z/(r2*xyDist);

        //Derivatives with respect to y.
        J[3]=y/r;
        J[4]=-x/x2y2Sum;
        J[5]=-y*z/(r2*xyDist);

        //Derivatives with respect to z.
        J[6]=z/r;
        J[7]=0;
        J[8]=xyDist/r2;
    } else {//Assume systemType==1
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
