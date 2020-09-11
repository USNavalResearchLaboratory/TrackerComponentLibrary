/*These function perform conversions from bistatic r-u-v coordinates into
 *Cartesian coordinates in different ways. See the Matlab implementation of
 *ruv2Cart for more insight on the algorithms.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include<CoordFuncs.hpp>
//For sqrt
#include <cmath>

void ruv2CartGenCPP(double *retData,const double *z,const bool useHalfRange,const double *zTx,const double *zRx,const double *M, bool hasW) {
/*RUV2CARTGENCPP A C++ function to convert a point in range, and direction
 *           cosines, possibly including w to Cartesian coordinates.
 *
 *INPUTS: retData A pointer to an array of doubles with 3 elements to
 *                hold the result in [x,y,z] order.
 *              z The 3X1 or 4X1 (if hasW is true) point in r-u-v-(w)
 *                coordinates to be converted, ordered [range;u;v;w].
 *   useHalfRange A boolean value specifying whether the bistatic (round-
 *                trip) range value has been divided by two. 
 *            zTx The 3X1 [x;y;z] location vector of the transmitter in
 *                global Cartesian coordinates.
 *            zRx The 3X1 [x;y;z] location vector of the receiver in global
 *                Cartesian coordinates.
 *             M  A 3X3  rotation matrices to go from the alignment of the
 *                global coordinate system to that at the receiver. It is
 *                stored one columns after the other, consistent with how
 *                Matlab indexes matrices.
 *           hasW A boolean value indicating whether the measurement has
 *                the third direction vector component (after u and v) and
 *                is thus 4X1 instrad of 3X1.
 *
 *OUTPUTS: None. The results are placed in retData.
 *
 *See the comments to the Matlab function ruv2Cart for more information
 *on how this function works.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/

    double rB,u,v,r1;
    double uvMag2;
    double uVec[3];
    double zTxL[3];
    
    //Extract the components
    rB=z[0];
    u=z[1];
    v=z[2];
    
    //The bistatic range is used in the conversions below.
    if(useHalfRange){
       rB=2*rB; 
    }
    
    if(hasW) {
        uVec[0]=u;
        uVec[1]=v;
        uVec[2]=z[3];
    } else {
        //If the magnitude is too large, normalize it so that one does
        //not try to take the square root of a negative number.
        uvMag2=u*u+v*v;
        if(uvMag2>1) {
            double uvMag=sqrt(uvMag2);
            u=u/uvMag;
            v=v/uvMag;

        //uVec is a unit vector pointing in the direction of the target from
        //the receiver.
            uVec[0]=u;
            uVec[1]=v;
            uVec[2]=0;
        } else {
            uVec[0]=u;
            uVec[1]=v;
            uVec[2]=sqrt(1.0-uvMag2);
        }
    }

    //Compute the transmitter location in the receiver's local coordinate
    //system. zTxL=M*(zTx-zRx);
    {
    double diff0=zTx[0]-zRx[0];
    double diff1=zTx[1]-zRx[1];
    double diff2=zTx[2]-zRx[2];
    zTxL[0]=M[0]*diff0+M[3]*diff1+M[6]*diff2;
    zTxL[1]=M[1]*diff0+M[4]*diff1+M[7]*diff2;
    zTxL[2]=M[2]*diff0+M[5]*diff1+M[8]*diff2;
    }
    
    r1=(rB*rB-(zTxL[0]*zTxL[0]+zTxL[1]*zTxL[1]+zTxL[2]*zTxL[2]))/(2.0*(rB-(uVec[0]*zTxL[0]+uVec[1]*zTxL[1]+uVec[2]*zTxL[2])));
    
    //This deals with the case where a point with zero range in the
    //monostatic case is passed. In the bistatic case, a zero range is
    //invalid.
    if(rB==0) {
        r1=0;
    }

    //zL is the Cartesian location in the local coordinate system of
    //the receiver
    {
    double zL[3];
    zL[0]=r1*uVec[0];
    zL[1]=r1*uVec[1];
    zL[2]=r1*uVec[2];
        
    //Convert to global Cartesian coordinates. The transpose of a rotation
    //matrix is its inverse. retData=M'*zL+zRx;
    retData[0]=M[0]*zL[0]+M[1]*zL[1]+M[2]*zL[2]+zRx[0];
    retData[1]=M[3]*zL[0]+M[4]*zL[1]+M[5]*zL[2]+zRx[1];
    retData[2]=M[6]*zL[0]+M[7]*zL[1]+M[8]*zL[2]+zRx[2];
    }  
}

/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
