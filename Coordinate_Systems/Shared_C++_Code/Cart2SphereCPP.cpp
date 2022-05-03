/*These function perform conversions from Cartesian coordinates to
 *spherical coordinates. See the Matlab implementation of Cart2Sphere for
 *more insight on the algorithms.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
//For sqrt and trigonometric functions.
#include <math.h>

void Cart2SphereGenCPP(double *retData,const double *cartPoints,const size_t systemType,const bool useHalfRange,const double *zTx,const double *zRx,const double *M) {
/*CART2SPHEREGENCPP A C++ function to convert Cartesian points to bistatic
 *            range, azimuth and elevation.
 *
 *INPUTS: retData A pointer to an array of doubles with 3 elements to
 *                hold the result in [range;azimuth;elevation]. order.
 *     cartPoints A pointer to the 3X1 Cartesian points [x;y;z] to be
 *                converted.
 *     systemType An integer specifying the axis from which the angles are
 *                measured. Possible values are
 *                0 Azimuth is measured counterclockwise from the x-axis in
 *                  the x-y plane. Elevation is measured up from the x-y
 *                  plane (towards the z-axis). This is consistent with
 *                  common spherical coordinate systems for specifying
 *                  longitude (azimuth) and geocentric latitude
 *                  (elevation).
 *                1 Azimuth is measured counterclockwise from the z-axis in
 *                  the z-x plane. Elevation is measured up from the z-x
 *                  plane (towards the y-axis). This is consistent with
 *                  some spherical coordinate systems that use the z-axis
 *                  as the boresight direction of the radar.
 *                2 This is the same as 0 except instead of being given
 *                  elevation, one desires the angle away from the z-axis,
 *                  which is (pi/2-elevation).
 *                3 This is the same as 0 except azimuth is measured
 *                  clockwise from the y-axis in the x-y plane instead of
 *                  counterclockwise from the x-axis. This coordinate
 *                  system often arises when given "bearings" in a local
 *                  East-North-Up coordinate system, where the bearing
 *                  directions are measured East of North.
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
 *
 *OUTPUTS: None. The results are placed in retData.
 *
 *See the comments to the Matlab function Cart2Sphere for more information
 *on how this function works.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
 
    double zCL[3], zTxL[3];
    double diff0,diff1,diff2,r1,r2;

    //Compute the target location in the receiver's coordinate system.
    diff0=cartPoints[0]-zRx[0];
    diff1=cartPoints[1]-zRx[1];
    diff2=cartPoints[2]-zRx[2];
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
    
    r1=sqrt(zCL[0]*zCL[0]+zCL[1]*zCL[1]+zCL[2]*zCL[2]);//Receiver->target.
    diff0=zCL[0]-zTxL[0];
    diff1=zCL[1]-zTxL[1];
    diff2=zCL[2]-zTxL[2];
    r2=sqrt(diff0*diff0+diff1*diff1+diff2*diff2);
    retData[0]=r1+r2;//The bistatic range.
    
    if(systemType==0||systemType==2||systemType==3) {
        //The special case for two zeros deals with the fact that the
        //standard C++ library will normally throw a domain error.
        if(zCL[1]==0&&zCL[0]==0) {
            retData[1]=0;//Azimuth
        } else {
            retData[1]=atan2(zCL[1],zCL[0]);//Azimuth
        }
        //The atan2 formulation is numerically more accurate than the asin
        //asin formulation.
        retData[2]=atan2(zCL[2],hypot(zCL[0],zCL[1]));//=asin(zCL[2]/r1) Elevation
        
        if(systemType==2) {
            double pi=acos(-1.0);
            retData[2]=pi/2.0-retData[2];
        }else if(systemType==3) {
            double pi=acos(-1.0);
            retData[1]=pi/2.0-retData[1];
        }
    } else {//Assume systemType=1
        if(zCL[2]==0&&zCL[0]==0) {
            retData[1]=0;//Azimuth
        } else {
            retData[1]=atan2(zCL[0],zCL[2]);//Azimuth
        }
        retData[2]=atan2(zCL[1],hypot(zCL[2],zCL[0]));//=asin(zCL[1]/r1) Elevation
    }

    if(useHalfRange) {
        retData[0]=retData[0]/2;
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
