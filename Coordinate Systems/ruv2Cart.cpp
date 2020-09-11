/*RUV2CART Convert points in bistatic r-u-v (or r-u-v-w) coordinates into
*         Cartesian coordinates. If r-u-v coordinates are used, the target
*         is assumed to be in front of the receiver (local z coordinate is
*         positive). r-u-v coordinates consist of a bistatic range and
*         direction cosines at the receiver. The "direction cosines" u and
*         v are just the x and y coordinates of a unit vector from the
*         receiver to the target in the coordinate system at the receiver.
*         This basically assumes that the boresight direction of the
*         receiver is the z axis. Assuming the target is in front of the
*         receiver, the third unit vector coordinate is not needed.
*         However, r-u-v-w coordinates include the third component.
*         For monostatic coordinate conversion where the range is a one-way
*         range, set useHalfRange=true and zTx and zRx to the same value.
*
*INPUTS: z A 3XN matrix of vectors with elements [r;u;v], where r is the
*          bistatic range from the transmitter to the target to the
*          receiver, and u and v are direction cosines. Each u,v pair
*          should have a magnitude less than or equal to one. If the
*          magnitude is greater than one, then the pair is normalized
*          before conversion to avoid imaginary results. Alternatively, one
*          can pass a 4XN matrix of [r;u;v;w] vectors where [u;v;w] form a
*          full unit vector in the receiver's local 3D Cartesian
*          coordinates.
* useHalfRange A boolean value specifying whether the bistatic range value
*          should be divided by two. This normally comes up when operating
*          in monostatic mode, so that the range reported is a one-way
*          range. The default if this parameter is not provided, or an
*          empty matrix is passed, is false.
*      zTx The 3XN [x;y;z] location vectors of the transmitters in global
*          Cartesian coordinates. If this parameter is omitted or an empty
*          matrix is passed, then the transmitters are assumed to be at the
*          origin. If only a single vector is passed, then the transmitter
*          location is assumed the same for all of the target states being
*          converted. zTx can have more than 3 rows; additional rows are
*          ignored.
*      zRx The 3XN [x;y;z] location vectors of the receivers in Cartesian
*          coordinates.  If this parameter is omitted or an empty matrix
*          is passed, then the receivers are assumed to be at the origin.
*          If only a single vector is passed, then the receiver location
*          is assumed the same for all of the target states being
*          converted. zRx can have more than 3 rows; additional rows are
*          ignored.
*        M A 3X3XN hypermatrix of the rotation matrices to go from the
*          alignment of the global coordinate system to that at the
*          receiver. The z-axis of the local coordinate system of the
*          receiver is the pointing direction of the receiver. If omitted
*          or an empty matrix is passed, then it is assumed that the local
*          coordinate system is aligned with the global and M=eye(3) --the
*          identity matrix is used. If only a single 3X3 matrix is passed,
*          then is is assumed to be the same for all of the N conversions.
*
*OUTPUTS: zC The 3XN matrix of the converted points in [x;y;z] Cartesian
*            coordinates.
*
*Basic u and v direction cosines do not specify which side of the radar the
*target is on. That is, they do not specify the sign of z. This
*performs the conversion assuming that z is positive in the local
*coordinate system of the receiver. Details of the conversion are given in
*[1].
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*zC=ruv2Cart(z,useHalfRange,zTx,zRx,M)
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
    
    //This is 3 if there is no w component and 4 is there is a w components
    //to the measurements.
    size_t pointOffset;
    bool hasW;
    
    mxArray *retMat;
    double *retData;
    //These two could be declared const, but that would just require extra
    //typecasting, since the return value of mxGetDoubles for the inputs is
    //not const and would have to be typecase, or these would have to be
    //typecast.
    double zeroVector[3]={0,0,0};
    double identMat[9]={1,0,0,0,1,0,0,0,1};
    bool useHalfRange;

    if(nrhs<1||nrhs>5) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    pointOffset=mxGetM(prhs[0]);
    N=mxGetN(prhs[0]);

    if((pointOffset!=3&&pointOffset!=4)||N<1) {
       mexErrMsgTxt("The points have the wrong dimensionality.");
       return;
    }
    
    hasW=(pointOffset==4);
    
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

    retMat=mxCreateDoubleMatrix(3,N,mxREAL);
    retData=mxGetDoubles(retMat);
    
    //Convert each of the measurements
    for(i=0;i<N;i++) {
        ruv2CartGenCPP(retData,points,useHalfRange,zTx,zRx,M,hasW);
        retData=retData+3;
        points=points+pointOffset;
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
