/**SPHER2CART Convert points from spherical coordinates to Cartesian
*             coordinates. If no range is specified, then the angles are
*             converted into Cartesian unit vectors. An option allows for
*             the angles to be specified from different axes. Optionally, a
*             bistatic range can be used when considering bistatic
*             measurements in a local spherical coordinate system. 
*
*INPUTS: points One or more points given in terms of range, azimuth and
*           elevation, with the angles in radians, or in terms of just
*           azimuth and elevation if Cartesian unit vectors are desired. To
*           convert N points, points is a 3XN matrix with each column
*           having the format [range;azimuth; elevation] or it is a 2XN
*           matrix with each column having format [azimuth; elevation] if
*           unit vectors are desired. Note that many math texts use a polar
*           angle (pi/2-elevation) in place of elevation. A polar angle is
*           also known as a colatitude, an inclination angle, a zenith
*           angle, and a normal angle. systemType=2 supports a polar angle.
* systemType An optional parameter specifying the axis from which the
*           angles are measured in radians. Possible values are
*           0 (The default if omitted) Azimuth is measured 
*             counterclockwise from the x-axis in the x-y plane. Elevation
*             is measured up from the x-y plane (towards the z-axis). This
*             is consistent with common spherical coordinate systems for
*             specifying longitude (azimuth) and geocentric latitude
*             (elevation).
*           1 Azimuth is measured counterclockwise from the z-axis in the
*             z-x plane. Elevation is measured up from the z-x plane
*             (towards the y-axis). This is consistent with some spherical
*             coordinate systems that use the z axis as the boresight
*             direction of the radar.
*           2 This is the same as 0 except instead of being given
*             elevation, one is given the angle away from the z-axis, which
*             is (pi/2-elevation).
*           3 This is the same as 0 except azimuth is measured clockwise
*             from the y-axis in the x-y plane instead of counterclockwise
*             from the x-axis. This coordinate system often arises when
*             given "bearings" in a local East-North-Up coordinate system,
*             where the bearing directions are measured East of North.
* useHalfRange An optional boolean value specifying whether the bistatic
*           (round-trip) range value has been divided by two. This normally
*           comes up when operating in monostatic mode (the most common
*           type of spherical coordinate system), so that the range
*           reported is a one-way range (or just half a bistatic range).
*           The default if this parameter is not provided is false if zTx
*           is provided and true if it is omitted (monostatic). If no
*           range values are provided, an empty matrix can be passed.
*       zTx The 3XN [x;y;z] location vectors of the transmitters in global
*           Cartesian coordinates. If this parameter is omitted, then the
*           transmitters are assumed to be at the origin. If only a single
*           vector is passed, then the transmitter location is assumed the
*           same for all of the target states being converted. zTx can have
*           more than 3 rows; additional rows are ignored. If monostatic or
*           no range values are provided, an empty matrix can be passed.
*       zRx The 3XN [x;y;z] location vectors of the receivers in Cartesian
*           coordinates.  If this parameter is omitted, then the
*           receivers are assumed to be at the origin. If only a single
*           vector is passed, then the receiver location is assumed the
*           same for all of the target states being converted. zRx can have
*           more than 3 rows; additional rows are ignored. If monostatic or
*           no range values are provided, an empty matrix can be passed.
*         M A 3X3XN hypermatrix of the rotation matrices to go from the
*           alignment of the global coordinate system to that at the
*           receiver. If omitted, then it is assumed that the local
*           coordinate system is aligned with the global and M=eye(3) --the
*           identity matrix is used. If only a single 3X3 matrix is passed,
*           then it is assumed to be the same for all of the N conversions.
* flipNegRange An optional boolean parameter. If true, the absolute value
*           of all range components, if present, is taken before
*           conversion.
*
*OUTPUTS: cartPoints For N points, cartPoints is a 3XN matrix of the
*                    converted points with each column having the format
*                    [x;y;z]. If no range values were passed, then all of
*                    the vectors are unit direction vectors.
*
*The conversion from spherical to Cartesian coordinates is given in [1].
*However, when considering the bistatic scenario, concepts from the
*bistatic r-u-v to Cartesian conversion described in [2] are used.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*cartPoints=spher2Cart(points,systemType,useHalfRange,zTx,zRx,M,flipNegRange)
*
*REFERENCES:
*[1] R. L. Duncombe, "Computational techniques," in Explanatory Supplement
*    to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
*    Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013,
*    ch. 14.4.4.1.
*[2] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
*    bistatic measurements," IEEE Aerospace and Electronic Systems
*    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
*
*February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "CoordFuncs.hpp"
#include <cmath>
        
void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double *points, *zTx, *zRx, *M;
    int systemType;
    size_t i, N, numRows;
    //If multiple values are passed for zTx, zRx or M, then the three
    //offsets below are used to move to the next value when going through
    //the measurements. However, if onyl a single value is passed, but
    //multiple measurements are present, then these will stay 0, causing
    //the single value to be reused. 
    size_t zTxOffset=0;
    size_t zRxOffset=0;
    size_t MOffset=0;
    
    mxArray *retMat;
    double *retData;
    //These two could be declared const, but that would just require extra
    //typecasting, since the return value of mxGetDoubles for the inputs is
    //not const and would have to be typecase, or these would have to be
    //typecast.
    double zeroVector[3]={0,0,0};
    double identMat[9]={1,0,0,0,1,0,0,0,1};
    bool useHalfRange;
    bool flipNegRange=false;

    if(nrhs<1||nrhs>7) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    if(nrhs>6&&!mxIsEmpty(prhs[6])) {
        flipNegRange=getBoolFromMatlab(prhs[6]);
    }
    
    numRows=mxGetM(prhs[0]);
    N=mxGetN(prhs[0]);

    if((numRows!=3&&numRows!=2)||N<1) {
       mexErrMsgTxt("The points have the wrong dimensionality.");
       return;
    }
    
    checkRealDoubleArray(prhs[0]);
    points=mxGetDoubles(prhs[0]);
    
    if(nrhs<2||mxIsEmpty(prhs[1])) {
        systemType=0;
    } else {
        systemType=getIntFromMatlab(prhs[1]);
        
        if(systemType!=0&&systemType!=1&&systemType!=2&&systemType!=3) {
            mexErrMsgTxt("Invalid system type specified.");
            return;
        }
    }

    //If the measurement have range components.
    if(numRows==3) {
        if(nrhs<3||mxIsEmpty(prhs[2])) {
            if(nrhs<4||mxIsEmpty(prhs[3])) {//If zTx is omitted
                useHalfRange=true;
            } else {
                useHalfRange=false; 
            }
        } else {
            useHalfRange=getBoolFromMatlab(prhs[2]);
        }
    
        if(nrhs<4||mxIsEmpty(prhs[3])) {
            zTx=zeroVector;
        } else {
            size_t numVecs, numTxRows;

            checkRealDoubleArray(prhs[3]);
            numVecs=mxGetN(prhs[3]);
            numTxRows=mxGetM(prhs[3]);

            if(numTxRows<3||numVecs<1||(numVecs!=N&&numVecs!=1)) {
                mexErrMsgTxt("The transmitter locations have the wrong dimensionality.");
                return;
            }
            
            zTx=mxGetDoubles(prhs[3]);
            if(numVecs==N) {
                zTxOffset=numTxRows;
            }
        }
    
        if(nrhs<5||mxIsEmpty(prhs[4])) {
            //Allocate a vector of the correct size of doubles.
            zRx=zeroVector;
        } else {
            size_t numVecs, numRxRows;

            checkRealDoubleArray(prhs[4]);
            numVecs=mxGetN(prhs[4]);
            numRxRows=mxGetM(prhs[4]);

            if(numRxRows<3||numVecs<1||(numVecs!=N&&numVecs!=1)) {
                mexErrMsgTxt("The transmitter locations have the wrong dimensionality.");
                return;
            }
            
            zRx=mxGetDoubles(prhs[4]);
            
            if(numVecs==N) {
                zRxOffset=numRxRows;
            }
        }
    }

    if(nrhs<6||mxIsEmpty(prhs[5])) {
        M=identMat;
    } else {
        size_t numMats;
        size_t numMatDims;
        const size_t *matDims;

        checkRealDoubleHypermatrix(prhs[5]);
        numMatDims=mxGetNumberOfDimensions(prhs[5]);
        matDims=mxGetDimensions(prhs[5]);

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
        
        M=mxGetDoubles(prhs[5]);
        if(numMats==N) {
            MOffset=9;
        }
    }

    retMat=mxCreateDoubleMatrix(3,N,mxREAL);
    retData=mxGetDoubles(retMat);
    
    //Convert each of the measurements
    if(numRows==3) {//There is range.
        if(flipNegRange) {
            for(i=0;i<N;i++) {
                double curPoint[3];
                curPoint[0]=fabs(points[0]);
                curPoint[1]=points[1];
                curPoint[2]=points[2];

                spher2CartGenCPP(retData,curPoint,systemType,useHalfRange,zTx,zRx,M);
                retData=retData+3;
                points=points+3;
                zTx+=zTxOffset;
                zRx+=zRxOffset;
                M=M+MOffset;
            }
        } else {
            for(i=0;i<N;i++) {
                spher2CartGenCPP(retData,points,systemType,useHalfRange,zTx,zRx,M);
                retData=retData+3;
                points=points+3;
                zTx+=zTxOffset;
                zRx+=zRxOffset;
                M=M+MOffset;
            }
        }
    } else {//M==2 --there is no range
        for(i=0;i<N;i++) {
            spher2CartNoRangeCPP(retData,points,systemType,M);

            retData=retData+3;
            points=points+2;
            M=M+MOffset;
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
