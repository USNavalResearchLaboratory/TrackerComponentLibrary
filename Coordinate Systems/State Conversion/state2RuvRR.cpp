/*STATE2RUVRR Convert state vectors consisting of at least 3D position and
*             velocity in 3D space into local bistatic r-u-v coordinates
*             with non-relativistic range rate.
*
*INPUTS: xTar An xDimXN matrix of N target states consisting of 3D position
*             and velocity components in the order
*             xTar=[xPosition;xVelocity] and possible other components,
*             which will be ignored.
* useHalfRange A boolean value specifying whether the bistatic range value
*             should be divided by two. This normally comes up when
*             operating in monostatic mode, so that the range reported is a
*             one-way range. The default if this parameter is not provided
*             is false.
*         xTx An xTxDimXN matrix of the states of the transmitters
*             consisting of stacked 3D position and velocity components.
*             Other components will be ignored. If this parameter is
*             omitted, the transmitters are assumed to be stationary at the
*             origin. If only a single vector is passed, then the
*             transmitter state is assumed the same for all of the target
*             states being converted.
*         xRx An xRxDimXN matrix of the states of the receivers consisting
*             of stacked 3D position and velocity components. Other
*             components will be ignored. If this parameter is omitted, the
*             receivers are assumed to be stationary at the origin. If only
*             a single vector is passed, then the receiver state is assumed
*             the same for all of the target states being converted.
*           M A 3X3XN hypermatrix of the rotation matrices to go from the
*             alignment of the global coordinate system to that at the
*             receiver. The z-axis of the local coordinate system of the
*             receiver is the pointing direction of the receiver. If
*             omitted, then it is assumed that the local coordinate system
*             is aligned with the global and M=eye(3) --the identity matrix
*             is used. If only a single 3X3 matrix is passed, then is is
*             assumed to be the same for all of the N conversions.
*    includeW An optional boolean value indicating whether a third
*             direction cosine component should be included. The u and v
*             direction cosines are two parts of a 3D unit vector.
*             Generally, one might assume that the target is in front of
*             the sensor, so the third component would be positive and is
*             not needed. However, the third component can be included if
*             ambiguity exists. The default if this parameter is omitted or
*             an empty matrix is passed is false.
*
*OUTPUTS: z A 4XN (or 5Xn if includeW is true) matrix of the target states
*           in xTar converted into bistatic r-u-v(-w)  and bistatic range
*           rate coordinates. If useHalfRange=true, then the r component is
*           half the bistatic range and the range rate is correspondingly
*           halved.
*
*Details of the conversions are given in [1].
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
* z=state2RuvRR(xTar,useHalfRange,xTx,xRx,M,includeW)
*
*REFERENCES:
*[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
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
    double *points, *xTx, *xRx, *M;
    size_t i, N;
    //If multiple values are passed for xTx, xRx or M, then the three
    //offsets below are used to move to the next value when going through
    //the measurements. However, if only a single value is passed, but
    //multiple measurements are present, then these will stay 0, causing
    //the single value to be reused. 
    size_t xTxOffset=0;
    size_t xRxOffset=0;
    size_t MOffset=0;
    //The size of the state. We skip over components that are other than
    //just position and velocity.
    size_t posOffset;

    mxArray *retMat;
    double *retData;
    size_t retVecOffsetPos;
    //These two could be declared const, but that would just require extra
    //typecasting, since the return value of mxGetDoubles for the inputs is
    //not const and would have to be typecase, or these would have to be
    //typecast.
    double zeroVector[6]={0,0,0,0,0,0};
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

    if(posOffset<6||N<1) {
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
        xTx=zeroVector;
    } else {
        size_t numVecs, numRows;

        checkRealDoubleArray(prhs[2]);
        numVecs=mxGetN(prhs[2]);
        numRows=mxGetM(prhs[2]);

        if(numRows<6||numVecs<1||(numVecs!=N&&numVecs!=1)) {
            mexErrMsgTxt("The transmitter locations have the wrong dimensionality.");
            return;
        }

        xTx=mxGetDoubles(prhs[2]);
        if(numVecs==N) {
            xTxOffset=numRows;
        }
    }

    if(nrhs<4||mxIsEmpty(prhs[3])) {
        //Allocate a vector of the correct size of doubles.
        xRx=zeroVector;
    } else {
        size_t numVecs, numRows;

        checkRealDoubleArray(prhs[3]);
        numVecs=mxGetN(prhs[3]);
        numRows=mxGetM(prhs[3]);

        if(numRows<6||numVecs<1||(numVecs!=N&&numVecs!=1)) {
            mexErrMsgTxt("The transmitter locations have the wrong dimensionality.");
            return;
        }

        xRx=mxGetDoubles(prhs[3]);

        if(numVecs==N) {
            xRxOffset=numRows;
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
        retVecOffsetPos=4;
        retMat=mxCreateDoubleMatrix(5,N,mxREAL);
    } else {
        retVecOffsetPos=3;
        retMat=mxCreateDoubleMatrix(4,N,mxREAL);
    }

    retData=mxGetDoubles(retMat);
    
    //Convert each of the measurements
    for(i=0;i<N;i++) {
        Cart2RuvGenCPP(retData,points,useHalfRange,xTx,xRx,M,includeW);
        retData=retData+retVecOffsetPos;
        *retData=getRangeRate3DCPP(points,useHalfRange,xTx,xRx);
        retData=retData+1;
        points=points+posOffset;
        xTx=xTx+xTxOffset;
        xRx=xRx+xRxOffset;
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
