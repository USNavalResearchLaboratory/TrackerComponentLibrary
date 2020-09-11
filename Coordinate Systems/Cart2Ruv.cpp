/*CART2RUV Convert points in Cartesian coordinates (either the global
%         system or a local system at the receiver) into local bistatic
%         r-u-v coordinates of the receiver. r-u-v coordinates consist of a
%         bistatic range and direction cosines at the receiver. The
%         "direction cosines" u and v are just the x and y coordinates of a
%         unit vector from the receiver to the target in the coordinate
%         system at the receiver. This basically assumes that the boresight
%         direction of the receiver is the z axis. Assuming the target is
%         in front of the receiver, the third unit vector coordinate is not
%         needed. However, with the includeW option, it can be provided,
%         resulting in r-u-v-w coordinates.
%
%INPUT: zC A 3XN matrix of Cartesian points in global [x;y;z] Cartesian
%          coordinates. Extra rows are ignored.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided (or an
%          empty matrix is provided) is false.
%      zTx The 3XN [x;y;z] location vectors of the transmitters in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the transmitters are assumed to be at the
%          origin. If only a single vector is passed, then the transmitter
%          location is assumed the same for all of the target states being
%          converted. zTx can have more than 3 rows; additional rows are
%          ignored.
%      zRx The 3XN [x;y;z] location vectors of the receivers in Cartesian
%          coordinates.  If this parameter is omitted or an empty matrix
%          is passed, then the receivers are assumed to be at the origin.
%          If only a single vector is passed, then the receiver location
%          is assumed the same for all of the target states being converted
%          zRx can have more than 3 rows; additional rows are ignored.
%        M A 3X3XN hypermatrix of the rotation matrices to go from the
%          alignment of the global coordinate system to that at the
%          receiver. The z-axis of the local coordinate system of the
%          receiver is the pointing direction of the receiver. If omitted
%          or an empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(3) --the
%          identity matrix is used. If only a single 3X3 matrix is passed,
%          then is is assumed to be the same for all of the N conversions.
% includeW An optional boolean value indicating whether a third direction
%          cosine component should be included. The u and v direction
%          cosines are two parts of a 3D unit vector. Generally, one might
%          assume that the target is in front of the sensor, so the third
%          component would be positive and is not needed. However, the
%          third component can be included if ambiguity exists. The default
%          if this parameter is omitted or an empty matrix is passed is
%          false.
%
%OUTPUT: z The 3XN (or 4XN if includeW is true) matrix of location vectors
%          of the points in bistatic [r;u;v] coordinates. If
%          useHalfRange=true, then the r component is half the bistatic
%          range (half the round-trip range for a monostatic scenario).
%
%Details of the conversion are given in [1].
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
* z=Cart2Ruv(zC,useHalfRange,zTx,zRx,M,includeW)
*
*REFERENCES:
*[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
*    bistatic measurements," IEEE Aerospace and Electronic Systems 
*    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
*
*April 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "CoordFuncs.hpp"
        
void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double *points, *zTx, *zRx, *M;
    size_t i, N;
    //If multiple values are passed for zTx, zRx or M, then the three
    //offsets below are used to move to the next value when going through
    //the measurements. However, if only a single value is passed, but
    //multiple measurements are present, then these will stay 0, causing
    //the single value to be reused. 
    size_t zTxOffset=0;
    size_t zRxOffset=0;
    size_t MOffset=0;
    //The size of x. We skip over components that are other than 3D
    //position.
    size_t posOffset;

    mxArray *retMat;
    double *retData;
    size_t retVecOffset;
    //These two could be declared const, but that would just require extra
    //typecasting, since the return value of mxGetDoubles for the inputs is
    //not const and would have to be typecase, or these would have to be
    //typecast.
    double zeroVector[3]={0,0,0};
    double identMat[9]={1,0,0,0,1,0,0,0,1};
    bool useHalfRange;
    bool includeW;

    if(nrhs<1||nrhs>6) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    N=mxGetN(prhs[0]);
    posOffset=mxGetM(prhs[0]);

    if(posOffset<3||N<1) {
       mexErrMsgTxt("The points have the wrong dimensionality.");
       return;
    }
    
    checkRealDoubleArray(prhs[0]);
    points=mxGetDoubles(prhs[0]);
    
    if(nrhs<2||mxIsEmpty(prhs[1])) {
        useHalfRange=true;
    } else {
        useHalfRange=getBoolFromMatlab(prhs[1]);
    }

    if(nrhs<3||mxIsEmpty(prhs[2])) {
        zTx=zeroVector;
    } else {
        size_t numVecs, numRows;

        checkRealDoubleArray(prhs[2]);
        numVecs=mxGetN(prhs[2]);
        numRows=mxGetM(prhs[2]);

        if(numRows<3||numVecs<1||(numVecs!=N&&numVecs!=1)) {
            mexErrMsgTxt("The transmitter locations have the wrong dimensionality.");
            return;
        }

        zTx=mxGetDoubles(prhs[2]);
        if(numVecs==N) {
            zTxOffset=numRows;
        }
    }

    if(nrhs<4||mxIsEmpty(prhs[3])) {
        //Allocate a vector of the correct size of doubles.
        zRx=zeroVector;
    } else {
        size_t numVecs, numRows;

        checkRealDoubleArray(prhs[3]);
        numVecs=mxGetN(prhs[3]);
        numRows=mxGetM(prhs[3]);

        if(numRows<3||numVecs<1||(numVecs!=N&&numVecs!=1)) {
            mexErrMsgTxt("The transmitter locations have the wrong dimensionality.");
            return;
        }

        zRx=mxGetDoubles(prhs[3]);

        if(numVecs==N) {
            zRxOffset=numRows;
        }
    }

    if(nrhs<5||mxIsEmpty(prhs[4])) {
        M=identMat;
    } else {
        size_t numMats;
        size_t numMatDims;
        const size_t *matDims;

        checkRealDoubleHypermatrix(prhs[4]);
        numMatDims=mxGetNumberOfDimensions(prhs[4]);
        matDims=mxGetDimensions(prhs[4]);

        if(numMatDims<2||numMatDims>3||matDims[0]!=3||matDims[1]!=3) {
           mexErrMsgTxt("The rotation matrices have the wrong dimensionality."); 
        }
        
        if(numMatDims==2) {
            numMats=1;
        } else {
            numMats=matDims[2];
        }
        
        if(numMats!=1&&numMats!=N) {
            mexErrMsgTxt("The rotation matrices have the wrong dimensionality."); 
        }
        
        M=mxGetDoubles(prhs[4]);
        if(numMats==N) {
            MOffset=9;
        }
    }

    if(nrhs<6||mxIsEmpty(prhs[5])) {
        includeW=false;
    } else {
        includeW=getBoolFromMatlab(prhs[5]);
    }
    
    if(includeW) {
        retVecOffset=4;
        retMat=mxCreateDoubleMatrix(4,N,mxREAL);
    } else {
        retVecOffset=3;
        retMat=mxCreateDoubleMatrix(3,N,mxREAL);
    }

    retData=mxGetDoubles(retMat);
    
    //Convert each of the measurements
    for(i=0;i<N;i++) {
        Cart2RuvGenCPP(retData,points,useHalfRange,zTx,zRx,M,includeW);
        retData=retData+retVecOffset;
        points=points+posOffset;
        zTx=zTx+zTxOffset;
        zRx=zRx+zRxOffset;
        M=M+MOffset;
    }

    plhs[0]=retMat;
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
