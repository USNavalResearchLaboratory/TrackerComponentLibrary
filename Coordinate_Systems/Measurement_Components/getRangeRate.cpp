/*GETRANGERATE Obtain the bistatic range rates of targets in the absence
*             of refraction under non-relativistic mechanics, ignoring
*             atmospheric effects, when the transmitter, target and
*             receiver are all moving. See the comments to the Matlab
*             native implementation for more details.
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

    mxArray *retMat;
    double *retData;
    //These two could be declared const, but that would just require extra
    //typecasting, since the return value of mxGetDoubles for the inputs is
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
    const size_t xDim=mxGetM(prhs[0]);

    if((xDim!=6&&xDim!=4&&xDim!=2)||N<1) {
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

        if(numRows<xDim||numVecs<1||(numVecs!=N&&numVecs!=1)) {
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

        if(numRows<xDim||numVecs<1||(numVecs!=N&&numVecs!=1)) {
            mexErrMsgTxt("The transmitter locations have the wrong dimensionality.");
            return;
        }

        zRx=mxGetDoubles(prhs[3]);

        if(numVecs==N) {
            zRxOffset=numRows;
        }
    }

    retMat=mxCreateDoubleMatrix(1,N,mxREAL);
    retData=mxGetDoubles(retMat);
    
    //Convert each of the measurements
    if(xDim==4) {//If it is in 2D
        for(i=0;i<N;i++) {
            *retData=getRangeRate2DCPP(points,useHalfRange,zTx,zRx);
            retData=retData+1;
            points=points+4;
            zTx=zTx+zTxOffset;
            zRx=zRx+zRxOffset;
        }
    } else if(xDim==6) {//It is in 3D
        for(i=0;i<N;i++) {
            *retData=getRangeRate3DCPP(points,useHalfRange,zTx,zRx);
            retData=retData+1;
            points=points+6;
            zTx=zTx+zTxOffset;
            zRx=zRx+zRxOffset;
        }
    } else{//xDim=2 -It is 1D.
        for(i=0;i<N;i++) {
            *retData=getRangeRate1DCPP(points,useHalfRange,zTx,zRx);
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
