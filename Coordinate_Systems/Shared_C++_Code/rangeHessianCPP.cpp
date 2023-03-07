/**RANGEHESSIANCPP A C++-only implementations of a function for computing
 *the Hessian of bistatic range.  See the Matlab equivalent for more
 *comments.
 *
 *June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

//Needed for sqrt
#include <cmath>
#include "CoordFuncs.hpp"

void rangeHessianGenCPP(const size_t numDim,double *H,const double *x,const bool useHalfRange,const double *lTx,const double *lRx) {
//numDim can be 1,2, or 3.

    if(numDim==1) {
        H[0]=0;
    } else if(numDim==2) {
        const double deltaRxx=x[0]-lRx[0];
        const double deltaRxx2=deltaRxx*deltaRxx;
        const double deltaRxy=x[1]-lRx[1];
        const double deltaRxy2=deltaRxy*deltaRxy;
        const double normDeltRx=sqrt(deltaRxx2+deltaRxy2);
        const double normDeltRx3=normDeltRx*normDeltRx*normDeltRx;
        const double deltaTxx=x[0]-lTx[0];
        const double deltaTxx2=deltaTxx*deltaTxx;
        const double deltaTxy=x[1]-lTx[1];
        const double deltaTxy2=deltaTxy*deltaTxy;
        const double normDeltTx=sqrt(deltaTxx2+deltaTxy2);
        const double normDeltTx3=normDeltTx*normDeltTx*normDeltTx;
        
        H[0]=-deltaRxx2/normDeltRx3+1/normDeltRx;
        H[1]=-deltaRxx*deltaRxy/normDeltRx3;
        H[3]=-deltaRxy2/normDeltRx3+1/normDeltRx;
        
        if(normDeltTx!=0) {
            //When the transmitter is the target.
            H[0]+= -deltaTxx2/normDeltTx3+1/normDeltTx;
            H[1]+= -deltaTxx*deltaTxy/normDeltTx3;
            H[3]+= -deltaTxy2/normDeltTx3+1/normDeltTx;
        }
        H[2]=H[1];

        if(useHalfRange) {
            size_t i;
            
            for(i=0;i<4;i++) {
                H[i]=H[i]/2;
            }
        }
    } else {//numDim==3
        const double deltaRxx=x[0]-lRx[0];
        const double deltaRxx2=deltaRxx*deltaRxx;
        const double deltaRxy=x[1]-lRx[1];
        const double deltaRxy2=deltaRxy*deltaRxy;
        const double deltaRxz=x[2]-lRx[2];
        const double deltaRxz2=deltaRxz*deltaRxz;
        const double normDeltRx=sqrt(deltaRxx2+deltaRxy2+deltaRxz2);
        const double normDeltRx3=normDeltRx*normDeltRx*normDeltRx;
        const double deltaTxx=x[0]-lTx[0];
        const double deltaTxy=x[1]-lTx[1];
        const double deltaTxz=x[2]-lTx[2];
        const double deltaTxx2=deltaTxx*deltaTxx;
        const double deltaTxy2=deltaTxy*deltaTxy;
        const double deltaTxz2=deltaTxz*deltaTxz;
        const double normDeltTx=sqrt(deltaTxx2+deltaTxy2+deltaTxz2);
        const double normDeltTx3=normDeltTx*normDeltTx*normDeltTx;
    
        //H(1,1) in Matlab
        H[0]=-deltaRxx2/normDeltRx3+1/normDeltRx;
        //H(2,1) in Matlab
        H[1]=-deltaRxx*deltaRxy/normDeltRx3;
        //H(3,1) in Matlab
        H[2]=-deltaRxx*deltaRxz/normDeltRx3;
        //H(2,2)
        H[4]=-deltaRxy2/normDeltRx3+1/normDeltRx;
        //H(3,2)
        H[5]=-deltaRxy*deltaRxz/normDeltRx3;
        //H(3,3)
        H[8]=-deltaRxz2/normDeltRx3+1/normDeltRx;
        
        if(normDeltTx!=0) {
            H[0]+= -deltaTxx2/normDeltTx3+1/normDeltTx;
            H[1]+= -deltaTxx*deltaTxy/normDeltTx3;
            H[2]+= -deltaTxx*deltaTxz/normDeltTx3;
            H[4]+= -deltaTxy2/normDeltTx3+1/normDeltTx;
            H[5]+= -deltaTxy*deltaTxz/normDeltTx3;
            H[8]+= -deltaTxz2/normDeltTx3+1/normDeltTx;
        }
        //H(1,2)
        H[3]=H[1];
        //H(1,3)
        H[6]=H[2];
        //H(2,3)
        H[7]=H[5];

        if(useHalfRange) {
            size_t i;
            
            for(i=0;i<9;i++) {
                H[i]=H[i]/2;
            }
        }
    }
}

void rangeHessianCPP(const size_t numDim,double *H,const double *x,const bool useHalfRange) {
//numDim can be 1,2, or 3. This function is monostatic-only. In this
//instance, the target and the transmitter are not collocated.

    if(numDim==1) {
        H[0]=0;
    } else if(numDim==2) {
        const double xC=x[0];
        const double xC2=xC*xC;
        const double yC=x[1];
        const double yC2=yC*yC;
        const double r=sqrt(xC2+yC2);
        const double r3=r*r*r;
        
        H[0]=-xC2/r3+1/r;
        H[1]=-xC*yC/r3;
        H[2]=H[1];
        H[3]=-yC2/r3+1/r;
        
        if(useHalfRange==false) {
            size_t i;
            
            for(i=0;i<4;i++) {
                H[i]=H[i]*2;
            }
        }
    } else {//numDim==3
        const double xC=x[0];
        const double xC2=xC*xC;
        const double yC=x[1];
        const double yC2=yC*yC;
        const double zC=x[2];
        const double zC2=zC*zC;
        const double r=sqrt(xC2+yC2+zC2);
        const double r3=r*r*r;
    
        //H(1,1) in Matlab
        H[0]=-xC2/r3+1/r;
        //H(2,1) in Matlab
        H[1]=-xC*yC/r3;
        //H(3,1) in Matlab
        H[2]=-xC*zC/r3;
        //H(1,2)
        H[3]=H[1];
        //H(2,2)
        H[4]=-yC2/r3+1/r;
        //H(3,2)
        H[5]=-yC*zC/r3;
        //H(1,3)
        H[6]=H[2];
        //H(2,3)
        H[7]=H[5];
        //H(3,3)
        H[8]=-zC2/r3+1/r;
        
        if(useHalfRange==false) {
            size_t i;
            
            for(i=0;i<9;i++) {
                H[i]=H[i]*2;
            }
        }
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
