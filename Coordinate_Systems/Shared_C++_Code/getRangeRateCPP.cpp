/*These function obtain a non-relativistic range rate value (ignoring
 *atmospheric propagation effects) from a Cartesian state consisting of
 *position and velocity components. See the Matlab implementation of
 *Cart2ruv for more insight on the algorithms.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include<CoordFuncs.hpp>
//For sqrt
#include <cmath>

double getRangeRate2DCPP(const double *xTar,bool useHalfRange,const double *xTx,const double *xRx) {
/*GETRANGERATE2DGENCPP A C++ function to convert a Cartesian state in 2D
 *        into a non-relativistic range rate, ignoring atmospheric effects.
 *
 *INPUTS: xTar The 4X1 Cartesian position and velocity vectors
 *             [x;y;xDot;yDot].
 *  useHalfRange A boolean value specifying whether the bistatic (round-
 *             trip) range value (and hence the range rate) has been
 *             divided by two. 
 *         xTx The 4X1 [x;y;xDot;yDot] position and velocity vector of
 *             the transmitter in global Cartesian coordinates.
 *         xRx The 4X1 [x;y;xDot;yDot] position and velocity vector of
 *             the receiver in global Cartesian coordinates.
 *
 *OUTPUTS: rr The range rate as a double.
 *
 *See the comments to the Matlab function getRangeRate for more information
 *on how this function works.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    
    double vt[2],vr[2],vi[2],dtr[2],dtl[2];
    double magVal,rr;

    vt[0]=xTar[2];
    vt[1]=xTar[3];

    vr[0]=xRx[2];
    vr[1]=xRx[3];
    
    vi[0]=xTx[2];
    vi[1]=xTx[3];
    
    dtr[0]=xTar[0]-xRx[0];
    dtr[1]=xTar[1]-xRx[1];
    //Normalize
    magVal=sqrt(dtr[0]*dtr[0]+dtr[1]*dtr[1]);
    dtr[0]/=magVal;
    dtr[1]/=magVal;
    
    dtl[0]=xTar[0]-xTx[0];
    dtl[1]=xTar[1]-xTx[1];
    //Normalize
    magVal=sqrt(dtl[0]*dtl[0]+dtl[1]*dtl[1]);

    if(magVal==0) {
        dtl[0]=0;
        dtl[1]=0;
    } else {
        dtl[0]/=magVal;
        dtl[1]/=magVal;
    }
    
    rr=(dtr[0]+dtl[0])*vt[0]+(dtr[1]+dtl[1])*vt[1]
            -dtr[0]*vr[0]-dtr[1]*vr[1]
            -dtl[0]*vi[0]-dtl[1]*vi[1];
    
    if(useHalfRange) {
        rr=rr/2;
    }
    return rr;
}

double getRangeRate3DCPP(const double *xTar,bool useHalfRange,const double *xTx,const double *xRx) {
/*GETRANGERATE3DGENCPP A C++ function to convert a Cartesian state in 3D
 *        into a non-relativistic range rate, ignoring atmospheric effects.
 *
 *INPUTS: xTar The 6X1 Cartesian position and velocity vectors
 *             [x;y;z;xDot;yDot;zDot].
 * useHalfRange A boolean value specifying whether the bistatic (round-
 *             trip) range value (and hence the range rate) has been
 *             divided by two. 
 *         xTx The 6X1 [x;y;z;xDot;yDot;zDot] position and velocity
 *             vector of the transmitter in global Cartesian coordinates.
 *         xRx The 6X1 [x;y;z;xDot;yDot;zDot] position and velocity
 *             vector of the receiver in global Cartesian coordinates.
 *
 *OUTPUTS: rr The range rate as a double.
 *
 *See the comments to the Matlab function getRangeRate for more information
 *on how this function works.
 *
 *April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    
    double vt[3],vr[3],vi[3],dtr[3],dtl[3];
    double magVal, rr;

    vt[0]=xTar[3];
    vt[1]=xTar[4];
    vt[2]=xTar[5];

    vr[0]=xRx[3];
    vr[1]=xRx[4];
    vr[2]=xRx[5];
    
    vi[0]=xTx[3];
    vi[1]=xTx[4];
    vi[2]=xTx[5];
    
    dtr[0]=xTar[0]-xRx[0];
    dtr[1]=xTar[1]-xRx[1];
    dtr[2]=xTar[2]-xRx[2];
    //Normalize
    magVal=sqrt(dtr[0]*dtr[0]+dtr[1]*dtr[1]+dtr[2]*dtr[2]);
    dtr[0]/=magVal;
    dtr[1]/=magVal;
    dtr[2]/=magVal;
    
    dtl[0]=xTar[0]-xTx[0];
    dtl[1]=xTar[1]-xTx[1];
    dtl[2]=xTar[2]-xTx[2];
    //Normalize
    magVal=sqrt(dtl[0]*dtl[0]+dtl[1]*dtl[1]+dtl[2]*dtl[2]);
    if(magVal==0) {
        dtl[0]=0;
        dtl[1]=0;
        dtl[2]=0;
    } else {
        dtl[0]/=magVal;
        dtl[1]/=magVal;
        dtl[2]/=magVal;
    }
    
    rr= (dtr[0]+dtl[0])*vt[0]+(dtr[1]+dtl[1])*vt[1]+(dtr[2]+dtl[2])*vt[2]
            -dtr[0]*vr[0]-dtr[1]*vr[1]-dtr[2]*vr[2]
            -dtl[0]*vi[0]-dtl[1]*vi[1]-dtl[2]*vi[2];
    
    if(useHalfRange) {
        rr=rr/2;
    }
    return rr;
}

double getRangeRate1DCPP(const double *xTar,bool useHalfRange,const double *xTx,const double *xRx) {
/*GETRANGERATE3DGENCPP A C++ function to convert a Cartesian state in 1D
 *        into a non-relativistic range rate, ignoring atmospheric effects.
 *
 *INPUTS: xTar The 2X1 Cartesian position and velocity vectors
 *             [x;xDot].
 * useHalfRange A boolean value specifying whether the bistatic (round-
 *             trip) range value (and hence the range rate) has been
 *             divided by two. 
 *         xTx The 2X1 [x;xDot] position and velocity vector of the
 *             transmitter in global Cartesian coordinates.
 *         xTx The 2X1 [x;xDot] position and velocity vector of the
 *             receiver in global Cartesian coordinates.
 *
 *OUTPUTS: rr The range rate as a double.
 *
 *See the comments to the Matlab function getRangeRate for more information
 *on how this function works.
 *
 *March 20122 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
    const double dtr=xTar[0]-xRx[0];
    const double dtl=xTar[0]-xTx[0];

    const double dtrRat=dtr/sqrt(dtr*dtr);//Effectively a signum function.
    double dtlRat;
    
    if(dtl==0) {
        dtlRat=0;
    } else {
        dtlRat=dtl/sqrt(dtl*dtl);
    }

    double rr=(dtrRat+dtlRat)*xTar[1]-dtrRat*xRx[1]-dtlRat*xTx[1];

    if(useHalfRange) {
        return rr/2.0;
    } else {
        return rr;
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
