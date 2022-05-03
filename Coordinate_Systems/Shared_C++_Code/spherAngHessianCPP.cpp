/**SPHERANGHESSIANCPP A C++-only implementation of a function for
 *          computing the Hessian of azimuth and elevation in spherical
 *          coordinates.  See the Matlab equivalent for more comments.
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//Needed for sqrt
#include <cmath>
#include "CoordFuncs.hpp"

void spherAngHessianGenCPP(double *HTotal,const double *xG,const size_t systemType,const double *lRx,const double *M) {
    double xyz[3];
    double H[18];
    
    //xLocal=M*bsxfun(@minus,xG(1:3,:),lRx(1:3));
    xyz[0]=M[0]*(xG[0]-lRx[0])+M[3]*(xG[1]-lRx[1])+M[6]*(xG[2]-lRx[2]);
    xyz[1]=M[1]*(xG[0]-lRx[0])+M[4]*(xG[1]-lRx[1])+M[7]*(xG[2]-lRx[2]);
    xyz[2]=M[2]*(xG[0]-lRx[0])+M[5]*(xG[1]-lRx[1])+M[8]*(xG[2]-lRx[2]);
    
    spherAngHessianCPP(H,xyz,systemType);

    //Adjust for the rotated coordinate system.
    //HTotal(:,:,1)=M'*H(:,:,1)*M;
    //HTotal(:,:,2)=M'*H(:,:,2)*M;
    HTotal[0]=H[0]*M[0]*M[0]+M[1]*((H[1]+H[3])*M[0]+H[4]*M[1])+((H[2]+H[6])*M[0]+(H[5]+H[7])*M[1])*M[2]+H[8]*M[2]*M[2];
    HTotal[1]=H[0]*M[0]*M[3]+H[3]*M[1]*M[3]+H[6]*M[2]*M[3]+H[1]*M[0]*M[4]+H[4]*M[1]*M[4]+H[7]*M[2]*M[4]+H[2]*M[0]*M[5]+H[5]*M[1]*M[5]+H[8]*M[2]*M[5];
    HTotal[2]=H[0]*M[0]*M[6]+H[3]*M[1]*M[6]+H[6]*M[2]*M[6]+H[1]*M[0]*M[7]+H[4]*M[1]*M[7]+H[7]*M[2]*M[7]+H[2]*M[0]*M[8]+H[5]*M[1]*M[8]+H[8]*M[2]*M[8];
    HTotal[3]=H[0]*M[0]*M[3]+H[1]*M[1]*M[3]+H[2]*M[2]*M[3]+H[3]*M[0]*M[4]+H[4]*M[1]*M[4]+H[5]*M[2]*M[4]+H[6]*M[0]*M[5]+H[7]*M[1]*M[5]+H[8]*M[2]*M[5];
    HTotal[4]=H[0]*M[3]*M[3]+M[4]*((H[1]+H[3])*M[3]+H[4]*M[4])+((H[2]+H[6])*M[3]+(H[5]+H[7])*M[4])*M[5]+H[8]*M[5]*M[5];
    HTotal[5]=H[0]*M[3]*M[6]+H[3]*M[4]*M[6]+H[6]*M[5]*M[6]+H[1]*M[3]*M[7]+H[4]*M[4]*M[7]+H[7]*M[5]*M[7]+H[2]*M[3]*M[8]+H[5]*M[4]*M[8]+H[8]*M[5]*M[8];
    HTotal[6]=H[0]*M[0]*M[6]+H[1]*M[1]*M[6]+H[2]*M[2]*M[6]+H[3]*M[0]*M[7]+H[4]*M[1]*M[7]+H[5]*M[2]*M[7]+H[6]*M[0]*M[8]+H[7]*M[1]*M[8]+H[8]*M[2]*M[8];
    HTotal[7]=H[0]*M[3]*M[6]+H[1]*M[4]*M[6]+H[2]*M[5]*M[6]+H[3]*M[3]*M[7]+H[4]*M[4]*M[7]+H[5]*M[5]*M[7]+H[6]*M[3]*M[8]+H[7]*M[4]*M[8]+H[8]*M[5]*M[8];
    HTotal[8]=H[0]*M[6]*M[6]+M[7]*((H[1]+H[3])*M[6]+H[4]*M[7])+((H[2]+H[6])*M[6]+(H[5]+H[7])*M[7])*M[8]+H[8]*M[8]*M[8];
    HTotal[9]=H[9]*M[0]*M[0]+M[1]*((H[10]+H[12])*M[0]+H[13]*M[1])+((H[11]+H[15])*M[0]+(H[14]+H[16])*M[1])*M[2]+H[17]*M[2]*M[2];
    HTotal[10]=H[9]*M[0]*M[3]+H[12]*M[1]*M[3]+H[15]*M[2]*M[3]+H[10]*M[0]*M[4]+H[13]*M[1]*M[4]+H[16]*M[2]*M[4]+H[11]*M[0]*M[5]+H[14]*M[1]*M[5]+H[17]*M[2]*M[5];
    HTotal[11]=H[9]*M[0]*M[6]+H[12]*M[1]*M[6]+H[15]*M[2]*M[6]+H[10]*M[0]*M[7]+H[13]*M[1]*M[7]+H[16]*M[2]*M[7]+H[11]*M[0]*M[8]+H[14]*M[1]*M[8]+H[17]*M[2]*M[8];
    HTotal[12]=H[9]*M[0]*M[3]+H[10]*M[1]*M[3]+H[11]*M[2]*M[3]+H[12]*M[0]*M[4]+H[13]*M[1]*M[4]+H[14]*M[2]*M[4]+H[15]*M[0]*M[5]+H[16]*M[1]*M[5]+H[17]*M[2]*M[5];
    HTotal[13]=H[9]*M[3]*M[3]+M[4]*((H[10]+H[12])*M[3]+H[13]*M[4])+((H[11]+H[15])*M[3]+(H[14]+H[16])*M[4])*M[5]+H[17]*M[5]*M[5];
    HTotal[14]=H[9]*M[3]*M[6]+H[12]*M[4]*M[6]+H[15]*M[5]*M[6]+H[10]*M[3]*M[7]+H[13]*M[4]*M[7]+H[16]*M[5]*M[7]+H[11]*M[3]*M[8]+H[14]*M[4]*M[8]+H[17]*M[5]*M[8];
    HTotal[15]=H[9]*M[0]*M[6]+H[10]*M[1]*M[6]+H[11]*M[2]*M[6]+H[12]*M[0]*M[7]+H[13]*M[1]*M[7]+H[14]*M[2]*M[7]+H[15]*M[0]*M[8]+H[16]*M[1]*M[8]+H[17]*M[2]*M[8];
    HTotal[16]=H[9]*M[3]*M[6]+H[10]*M[4]*M[6]+H[11]*M[5]*M[6]+H[12]*M[3]*M[7]+H[13]*M[4]*M[7]+H[14]*M[5]*M[7]+H[15]*M[3]*M[8]+H[16]*M[4]*M[8]+H[17]*M[5]*M[8];
    HTotal[17]=H[9]*M[6]*M[6]+M[7]*((H[10]+H[12])*M[6]+H[13]*M[7])+((H[11]+H[15])*M[6]+(H[14]+H[16])*M[7])*M[8]+H[17]*M[8]*M[8];
}

void spherAngHessianCPP(double *H,const double *xG,const size_t systemType) {
    double x,y,z, x2, y2, z2, x4;
  
    x=xG[0];
    y=xG[1];
    z=xG[2];
    
    x2=x*x;
    y2=y*y;
    z2=z*z;
    x4=x2*x2;
    
    if(systemType==0) {
        double rxy, r4, r4xy, rxy3;
        
        r4xy=x2+y2;
        rxy=sqrt(r4xy);
        r4=(r4xy+z2);
        r4=r4*r4;
        r4xy=r4xy*r4xy;
        rxy3=rxy*rxy*rxy;

        //dAzdxdx
        H[0]=2*x*y/r4xy;
        //dAzdydy
        H[4]=-H[0];
        //dAzdzdz
        H[8]=0;
        //dAzdxdy
        H[3]=(y2-x2)/r4xy;
        //dAzdydx
        H[1]=H[3];
        //dAzdxdz
        H[6]=0;
        //dAzdzdx
        H[2]=H[6];
        //dAzdydz
        H[7]=0;
        //dAzdzdy
        H[5]=H[7];

        //dEldxdx
        H[9]=(z*(2*x4+x2*y2-y2*(y2+z2)))/(rxy3*r4);
        //dEldydy
        H[13]=-((z*(x4-2*y2*y2+x2*(z2-y2)))/(rxy3*r4));
        //dEldzdz
        H[17]=-(2*rxy*z)/r4;
        //dEldxdy
        H[12]=(x*y*z*(3*rxy*rxy+z2))/(rxy3*r4);
        //dEldydx
        H[10]=H[12];
        //dEldxdz
        H[15]=-((x*(x2+y2-z2))/(rxy*r4));
        //dEldzdx
        H[11]=H[15];
        //dEldydz
        H[16]=-((y*(x2+y2-z2))/(rxy*r4));
        //dEldzdy
        H[14]=H[16];
    } else if(systemType==1) {
        double x4,rxz, r4xz, rxz3, r4;

        x4=x2*x2;

        r4xz=x2+z2;
        rxz=sqrt(r4xz);
        rxz3=rxz*rxz*rxz;
        r4=(r4xz+y2);
        r4xz=r4xz*r4xz;
        r4=r4*r4;

        //dAzdxdx
        H[0]=-(2*x*z)/r4xz;
        //dAzdydy
        H[4]=0;
        //dAzdzdz
        H[8]=(2*x*z)/r4xz;
        //dAzdxdy
        H[3]=0;
        //dAzdydx
        H[1]=H[3];
        //dAzdxdz
        H[6]=(x2-z2)/r4xz;
        //dAzdzdx
        H[2]=H[6];
        //dAzdydz
        H[7]=0;
        //dAzdzdy
        H[5]=H[7];

        //dEldxdx
        H[9]=(y*(2*x4+x2*z2-z2*(y2+z2)))/(rxz3*r4);
        //dEldydy
        H[13]=-((2*y*rxz)/r4);
        //dEldzdz
        H[17]=-((y*(x4-2*z2*z2+x2*(y-z)*(y+z)))/(rxz3*r4));
        //dEldxdy
        H[12]=-((x*(x2-y2+z2))/(rxz*r4));
        //dEldydx
        H[10]=H[12];
        //dEldxdz
        H[15]=(x*y*z*(3*x2+y2+3*z2))/(rxz3*r4);
        //dEldzdx
        H[11]=H[15];
        //dEldydz
        H[16]=-((z*(x2-y2+z2))/(rxz*r4));
        //dEldzdy
        H[14]=H[16];
    } else if(systemType==3) {
        double rxy, r4, r4xy, rxy3;
        
        r4xy=x2+y2;
        rxy=sqrt(r4xy);
        r4=(r4xy+z2);
        r4=r4*r4;
        r4xy=r4xy*r4xy;
        rxy3=rxy*rxy*rxy;

        //dAzdxdx
        H[0]=-2*x*y/r4xy;
        //dAzdydy
        H[4]=-H[0];
        //dAzdzdz
        H[8]=0;
        //dAzdxdy
        H[3]=(x-y)*(x+y)/r4xy;
        //dAzdydz
        H[1]=H[3];
        //dAzdxdz
        H[6]=0;
        //dAzdzdx
        H[2]=H[6];
        //dAzdydz
        H[7]=0;
        //dAzdzdy
        H[5]=H[7];

        //dEldxdx
        H[9]=(z*(2*x4+x2*y2-y2*(y2+z2)))/(rxy3*r4);
        //dEldydy
        H[13]=-((z*(x4-2*y2*y2+x2*(z2-y2)))/(rxy3*r4));
        //dEldzdz
        H[17]=-(2*rxy*z)/r4;
        //dEldxdy
        H[12]=(x*y*z*(3*rxy*rxy+z2))/(rxy3*r4);
        //dEldydz
        H[10]=H[12];
        //dEldxdz
        H[15]=-((x*(x2+y2-z2))/(rxy*r4));
        //dEldzdx
        H[11]=H[15];
        //dEldydz
        H[16]=-((y*(x2+y2-z2))/(rxy*r4));
        //dEldzdy
        H[14]=H[16];
    }else {//systemType==2
        double r4xy, r4, rxy, x4, y4,rxy3;
        
        r4xy=x2+y2;
        r4=(r4xy+z2);
        r4=r4*r4;
        rxy=sqrt(r4xy);
        r4xy=r4xy*r4xy;
        rxy3=rxy*rxy*rxy;

        x4=x2*x2;
        y4=y2*y2;

        //dAzdxdx
        H[0]=(2*x*y)/r4xy;
        //dAzdydy
        H[4]=-H[0];
        //dAzdzdz
        H[8]=0;
        //dAzdxdy
        H[3]=(y2-x2)/r4xy;
        //dAzdydx
        H[1]=H[3];
        //dAzdxdz
        H[6]=0;
        //dAzdzdx
        H[2]=H[6];
        //dAzdydz
        H[7]=0;
        //dAzdzdy
        H[5]=H[7];

        //dEldxdx
        H[9]=(z*(-2*x4-x2*y2+y4+y2*z2))/(rxy3*r4);
        //dEldydy
        H[13]=(z*(x4-2*y4+x2*(z2-y2)))/(rxy3*r4);
        //dEldzdz
        H[17]=(2*rxy*z)/r4;
        //dEldxdy
        H[12]=-((x*y*z*(3*rxy*rxy+z2))/(rxy3*r4));
        //dEldydx
        H[10]=H[12];
        //dEldxdz
        H[15]=(x*(x2+y2-z2))/(rxy*r4);
        //dEldzdx
        H[11]=H[15];
        //dEldydz
        H[16]=(y*(x2+y2-z2))/(rxy*r4);
        //dEldzdy
        H[14]=H[16];
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
