/**SPHER2CARTCPP These function perform conversions from spherical
 *   coordinates into Cartesian coordinates in different ways. See the
 *   Matlab implementation of spher2Cart for more insight on the
 *   algorithms.
 *
 *February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
//For sin and cos
#include <math.h>

void spher2CartCPP(double *cartPoint,const double *point, size_t systemType) {
/*SPHER2CARTCPP A C++ function to convert a point from unrotated monostatic
 *           spherical coordinates with the sensor at the origin measing a
 *           one-way range to Cartesian coordinates
 *
 *INPUTS: cartPoint A pointer to an array of doubles with 3 elements to
 *                  hold the result in [x,y,z] order.
 *            point The 3X1 point in spherical coordinates to be converted,
 *                  ordered [range;azimuth;elevation].
 *       systemType An integer specifying the axis from which the angles
 *                  are measured. Possible values are
 *                   0 Azimuth is measured counterclockwise from the x-axis
 *                     in the x-y plane. Elevation is measured up from the
 *                     x-y plane (towards the z-axis). This is consistent
 *                     with common spherical coordinate systems for
 *                     specifying longitude (azimuth) and geocentric
 *                     latitude (elevation).
 *                   1 Azimuth is measured counterclockwise from the z-axis
 *                     in the z-x plane. Elevation is measured up from the
 *                     z-x plane (towards the y-axis). This is consistent
 *                     with some spherical coordinate systems that use the
 *                     z-axis as the boresight direction of the radar.
 *                   2 This is the same as 0 except instead of being given
 *                     elevation, one is given the angle away from the
 *                     z-axis, which is (pi/2-elevation).
 *                   3 This is the same as 0 except azimuth is measured
 *                     clockwise from the y-axis in the x-y plane instead
 *                     of counterclockwise from the x-axis. This coordinate
 *                     system often arises when given "bearings" in a local
 *                     East-North-Up coordinate system, where the bearing
 *                     directions are measured East of North.
 *
 *OUTPUTS: None. The results are placed in cartPoint.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    
    double r, azimuth, elevation;
    double cosEl;
    r=point[0];
    azimuth=point[1];
    elevation=point[2];
    
    if(systemType==2) {
        double pi=acos(-1.0);
        elevation=pi/2.0-elevation;
        systemType=0;
    } else if(systemType==3) {
        double pi=acos(-1.0);
        azimuth=pi/2.0-azimuth;
        systemType=0;
    }
    
    cosEl=cos(elevation);
    
    if(systemType==0) {
        cartPoint[0]=r*cos(azimuth)*cosEl;
        cartPoint[1]=r*sin(azimuth)*cosEl;
        cartPoint[2]=r*sin(elevation);
    } else {
        cartPoint[0]=r*sin(azimuth)*cosEl;
        cartPoint[1]=r*sin(elevation);
        cartPoint[2]=r*cos(azimuth)*cosEl;
    }
}

void spher2CartGenCPP(double *retData,const double *point,size_t systemType,const bool useHalfRange,const double *zTx,const double *zRx,const double *M) {
/*SPHER2CARTGENCPP A C++ function to convert a point in range, azimuth and
 *            elevation from bistatic spherical coordinates to Cartesian
 *            coordinates.
 *
 *INPUTS: retData A pointer to an array of doubles with 3 elements to
 *                hold the result in [x,y,z] order.
 *          point The 3X1 point in spherical coordinates to be converted,
 *                ordered [range;azimuth;elevation].
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
 *                  elevation, one is given the angle away from the z-axis,
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
 *See the comments to the Matlab function spher2Cart for more information
 *on how this function works.
 *
 *February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    double r, azimuth, elevation;
    double r1;
    double uVecL[3];
    double cosEl;
    double zTxL[3];
    double zL[3];
    
    r=point[0];
    
    if(useHalfRange==true) {
        r=2*r;
    }
    azimuth=point[1];
    elevation=point[2];
    
    if(systemType==2) {
        double pi=acos(-1.0);
        elevation=pi/2.0-elevation;
        systemType=0;
    } else if(systemType==3) {
        double pi=acos(-1.0);
        azimuth=pi/2.0-azimuth;
        systemType=0;
    }
    
    cosEl=cos(elevation);

    if(systemType==0) {
        uVecL[0]=cos(azimuth)*cosEl;
        uVecL[1]=sin(azimuth)*cosEl;
        uVecL[2]=sin(elevation);
    } else {//Assume systemType==1
        uVecL[0]=sin(azimuth)*cosEl;
        uVecL[1]=sin(elevation);
        uVecL[2]=cos(azimuth)*cosEl;
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
    
    //Compute the range from the transmitter to the target.
    //r1=(r^2-norm(zTxL)^2)/(2*(r-dot(uVec,zTxL)));
    {
    double num=r*r-(zTxL[0]*zTxL[0]+zTxL[1]*zTxL[1]+zTxL[2]*zTxL[2]);
    double denom=2*(r-uVecL[0]*zTxL[0]-uVecL[1]*zTxL[1]-uVecL[2]*zTxL[2]);
    
    r1=num/denom;
    }

    //Compute the Cartesian location in the local coordinate system of
    //the receiver.
    zL[0]=r1*uVecL[0];
    zL[1]=r1*uVecL[1];
    zL[2]=r1*uVecL[2];
    
    //Convert to global Cartesian coordinates. The transpose of a rotation
    //matrix is its inverse. retData=M'*zL+zRx;
    retData[0]=M[0]*zL[0]+M[1]*zL[1]+M[2]*zL[2]+zRx[0];
    retData[1]=M[3]*zL[0]+M[4]*zL[1]+M[5]*zL[2]+zRx[1];
    retData[2]=M[6]*zL[0]+M[7]*zL[1]+M[8]*zL[2]+zRx[2];
}

void spher2CartNoRangeCPP(double *retData,const double *point,size_t systemType,const double *M) {
/*SPHER2CARTNORANGECPP A C++ function to convert a direction in azimuth and
 *            angle from spherical coordinates to Cartesian coordinates.
 *
 *INPUTS: retData A pointer to an array of doubles with 3 elements to
 *                hold the result in [x,y,z] order.
 *          point The 2X1 point in spherical coordinates to be converted,
 *                ordered [azimuth;elevation].
 *     systemType An integer specifying the axis from which the angles are
 *                measured. Possible vaules are
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
 *                  elevation, one is given the angle away from the z-axis,
 *                  which is (pi/2-elevation).
 *                3 This is the same as 0 except azimuth is measured
 *                  clockwise from the y-axis in the x-y plane instead of
 *                  counterclockwise from the x-axis. This coordinate
 *                  system often arises when given "bearings" in a local
 *                  East-North-Up coordinate system, where the bearing
 *                  directions are measured East of North.
 *              M A 3X3  rotation matrices to go from the alignment of the
 *                global coordinate system to that at the receiver. It is
 *                stored one columns after the other, consistent with how
 *                Matlab indexes matrices.
 *
 *OUTPUTS: None. The results are placed in retData.
 *
 *See the comments to the Matlab function spher2Cart for more information
 *on how this function works.
 *
 *February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/  
    double azimuth, elevation;
    double cosEl;
    double uVecL[3];
    
    azimuth=point[0];
    elevation=point[1];
    
    if(systemType==2) {
        double pi=acos(-1.0);
        elevation=pi/2.0-elevation;
        systemType=0;
    } else if(systemType==3) {
        double pi=acos(-1.0);
        azimuth=pi/2.0-azimuth;
        systemType=0;
    }
    
    cosEl=cos(elevation);

    if(systemType==0) {
        uVecL[0]=cos(azimuth)*cosEl;
        uVecL[1]=sin(azimuth)*cosEl;
        uVecL[2]=sin(elevation);
    } else {
        uVecL[0]=sin(azimuth)*cosEl;
        uVecL[1]=sin(elevation);
        uVecL[2]=cos(azimuth)*cosEl;
    }
    
    //Convert to global coordinates.
    //retData=M'*uVecL;
    retData[0]=M[0]*uVecL[0]+M[1]*uVecL[1]+M[2]*uVecL[2];
    retData[1]=M[3]*uVecL[0]+M[4]*uVecL[1]+M[5]*uVecL[2];
    retData[2]=M[6]*uVecL[0]+M[7]*uVecL[1]+M[8]*uVecL[2];
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
