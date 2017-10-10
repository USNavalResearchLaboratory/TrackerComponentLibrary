/*GETRANGERATE Obtain the bistatic range rates of targets in the absence
*             of refraction under non-relativistic mechanics, ignoring
*             atmospheric effects, when the transmitter, target and
*             receiver are all moving.
*
*INPUTS: xTar The Cartesian state of the targets in 2D or 3D Cartesian
*             space. If 3D, then xTar is a 6XN matrix, (or where the first
*             3 components are position and the second 3 are velocity.
*             Otherwise, xTar is a 4XN matrix, where the first two
*             components are position and the second two velocity.
*             Components are ordered position (x,y,z) and velocity
*             (xDot,yDot,zDot).
* useHalfRange A boolean value specifying whether the bistatic range value
*             should be divided by two, which means that the range rate is
*             divided by two. This normally comes up when operating in
*             monostatic mode, so that the range reported is a one-way
*             range. The default if this parameter is not provided is
*             false.
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
*
*OUTPUTS: rr The 1XN bistatic range rates of the targets. If
*            useHalfRange=true, then the range rate is halved to reflect a
*            halved range.
*
*This assumes that the target state that is provided has the same
*dimensionality as the states of the transmitter and receiver and that all
*of the states are in Cartesian coordinates.
*
*A derivation of this non-relativistic approximation is given in [1].
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
* zrr=getRangeRate(xTar,useHalfRange,xTx,xRx)
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
    double *points, *zTx, *zRx;
    size_t i, N;
    //If multiple values are passed for zTx, zRx or M, then the three
    //offsets below are used to move to the next value when going through
    //the measurements. However, if only a single value is passed, but
    //multiple measurements are present, then these will stay 0, causing
    //the single value to be reused. 
    size_t zTxOffset=0;
    size_t zRxOffset=0;
    size_t xTarOffset;

    mxArray *retMat;
    double *retData;
    //These two could be declared const, but that would just require extra
    //typecasting, since the return value of mxGetData for the inputs is
    //not const and would have to be typecase, or these would have to be
    //typecast.
    double zeroVector[6]={0,0,0,0,0,0};
    bool useHalfRange;

    if(nrhs<1||nrhs>4) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    N=mxGetN(prhs[0]);
    xTarOffset=mxGetM(prhs[0]);

    if((mxGetM(prhs[0])!=6&&mxGetM(prhs[0])!=4)||N<1) {
       mexErrMsgTxt("The points have the wrong dimensionality.");
       return;
    }
    
    checkRealDoubleArray(prhs[0]);
    points=reinterpret_cast<double*>(mxGetData(prhs[0]));
    
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

        if(numRows<xTarOffset||numVecs<1||(numVecs!=N&&numVecs!=1)) {
            mexErrMsgTxt("The transmitter locations have the wrong dimensionality.");
            return;
        }

        zTx=reinterpret_cast<double*>(mxGetData(prhs[2]));
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

        if(numRows<xTarOffset||numVecs<1||(numVecs!=N&&numVecs!=1)) {
            mexErrMsgTxt("The transmitter locations have the wrong dimensionality.");
            return;
        }

        zRx=reinterpret_cast<double*>(mxGetData(prhs[3]));

        if(numVecs==N) {
            zRxOffset=numRows;
        }
    }

    retMat=mxCreateDoubleMatrix(1,N,mxREAL);
    retData=reinterpret_cast<double*>(mxGetData(retMat));
    
    //Convert each of the measurements
    if(xTarOffset==4) {//If it is in 2D
        for(i=0;i<N;i++) {
            *retData=getRangeRate2DCPP(points,useHalfRange,zTx,zRx);
            retData=retData+1;
            points=points+4;
            zTx=zTx+zTxOffset;
            zRx=zRx+zRxOffset;
        }
    } else {//It is in 3D
        for(i=0;i<N;i++) {
            *retData=getRangeRate3DCPP(points,useHalfRange,zTx,zRx);
            retData=retData+1;
            points=points+6;
            zTx=zTx+zTxOffset;
            zRx=zRx+zRxOffset;
        }
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
