/*These function perform conversions from Cartesian coordinates into
 *bistatic r-u-v coordinates. See the Matlab implementation of Cart2ruv for
 *more insight on the algorithms.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include<CoordFuncs.hpp>
//For sqrt
#include <cmath>

void Cart2RuvGenCPP(double *retData,const double *zC,bool useHalfRange,double *zTx,double *zRx,double *M,bool includeW) {
/*CART2RUVGENCPP A C++ function to convert a Cartesian point into range,
 *           and direction cosines, possibly including the w component.
 *
 *INPUTS: retData A pointer to an array of doubles with 3 elements to
 *                hold the result in [r;u;v] order or with 4 elements if
 *                includeW is true to hold [r;u;v;w].
 *             zC The 3X1 Cartesian points [x;y;z] to be converted.
 *   useHalfRange A boolean value specifying whether the bistatic (round-
 *                trip) range value has been divided by two. 
 *            zTx The 3X1 [x;y;z] location vector of the transmitter in
 *                global Cartesian coordinates.
 *            zRx The 3X1 [x;y;z] location vector of the receiver in global
 *                Cartesian coordinates.
 *             M  A 3X3 rotation matrices to go from the alignment of the
 *                global coordinate system to that at the receiver. It is
 *                stored one columns after the other, consistent with how
 *                Matlab indexes matrices.
 *       includeW A boolean value indicating whether retData has space for
 *                a fourth component and the fourth component should be
 *                included.
 *
 *OUTPUTS: None. The results are placed in retData.
 *
 *See the comments to the Matlab function Cart2ruv for more information
 *on how this function works.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/

    double zCL[3],zTxL[3];
    double r1,r2;
    double diff0,diff1,diff2;
 
    //Compute the target location in the receiver's coordinate system.
    diff0=zC[0]-zRx[0];
    diff1=zC[1]-zRx[1];
    diff2=zC[2]-zRx[2];
    zCL[0]=M[0]*diff0+M[3]*diff1+M[6]*diff2;
    zCL[1]=M[1]*diff0+M[4]*diff1+M[7]*diff2;
    zCL[2]=M[2]*diff0+M[5]*diff1+M[8]*diff2;

    //Compute the transmitter location in the receiver's local coordinate
    //system.
    diff0=zTx[0]-zRx[0];
    diff1=zTx[1]-zRx[1];
    diff2=zTx[2]-zRx[2];
    zTxL[0]=M[0]*diff0+M[3]*diff1+M[6]*diff2;
    zTxL[1]=M[1]*diff0+M[4]*diff1+M[7]*diff2;
    zTxL[2]=M[2]*diff0+M[5]*diff1+M[8]*diff2;

    r1=sqrt(zCL[0]*zCL[0]+zCL[1]*zCL[1]+zCL[2]*zCL[2]);//Receiver to target.

    diff0=zCL[0]-zTxL[0];
    diff1=zCL[1]-zTxL[1];
    diff2=zCL[2]-zTxL[2];
    r2=sqrt(diff0*diff0+diff1*diff1+diff2*diff2);//Target to transmitter.   

    retData[0]=r1+r2;
    retData[1]=zCL[0]/r1;
    retData[2]=zCL[1]/r1;

    if(includeW) {
        retData[3]=zCL[2]/r1;
    }

    if(useHalfRange) {
        retData[0]=retData[0]/2;
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
