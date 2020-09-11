/**GCRS2CIRS Convert vectors of position and possibly velocity from the
 *           Geocentric Celestrial Reference System (GCRS), a type of
 *           Earth-Centered Inertial (ECI) coordinate system, to the
 *           Celestial Intermediate Reference System (CIRS), which is
 *           offset by the rotation of the Celestial Intermediate Pole
 *           (CIP). The velocity conversion omits the centrifugal effects
 *           of the CIP motion, which have a period on the order of 14
 *           months and are thus small.
 *
 *INPUTS: x The NXnumVec collection of vectors to convert. N can be 3, or
 *          6. If the vectors are 3D, then are position. 6D vectors are
 *          assumed to be position and velocity.
 * Jul1, Jul2 Two parts of a Julian date given in terrestrial time (TT).
 *          The units of the date are days. The full date is the sum of
 *          both terms. The date is broken into two parts to provide more
 *          bits of precision. It does not matter how the date is split.
 *     dXdY dXdY=[dX;dY] are the celestial pole offsets with respect to
 *          the IAU 2006/2000A precession/nutation model in radians If this
 *          parameter is omitted, the value from the function getEOP will
 *          be used.
 *
 *OUTPUTS: vec A 3XN or 6XN matrix of vectors converted from GCRS
 *             coordinates to CIRS coordinates.
 *      rotMat The 3X3 rotation matrix used for the rotation of the
 *             positions and velocities.
 *
 *The conversion functions from the International Astronomical Union's
 *(IAU) Standard's of Fundamental Astronomy library are put together to get
 *the necessary rotation matrix for the position.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[vec,rotMat]=GCRS2CIRS(x,Jul1,Jul2);
 *or if more parameters are known,
 *[vec,rotMat]=GCRS2CIRS(x,Jul1,Jul2,dXdY);
 *
 *Different celestial coordinate systems are compared in [1].
 *
 *REFERENCES:
 *[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
 *    Temporal Coordinate Systems for Target Tracking," Formal Report,
 *    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
 *    173 pages.
 *
 *April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"
//For sqrt
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double TT1, TT2, dX, dY, *xVec;
    size_t numRow, numVec;
    mxArray *retMat;
    double *retData;
    double GCRS2CIRS[3][3];
    
    if(nrhs<3||nrhs>4){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    checkRealDoubleArray(prhs[0]);
    
    numRow = mxGetM(prhs[0]);
    numVec = mxGetN(prhs[0]);
    
    if(!(numRow==3||numRow==6)) {
        mexErrMsgTxt("The input vector has a bad dimensionality.");
    }

    xVec=mxGetDoubles(prhs[0]);
    TT1=getDoubleFromMatlab(prhs[1]);
    TT2=getDoubleFromMatlab(prhs[2]);
    
    //If some values from the function getEOP will be needed.
    if(nrhs<4||mxIsEmpty(prhs[3])) {
        mxArray *retVals[2];
        double *dXdY;
        mxArray *JulUTCMATLAB[2];
        double JulUTC[2];
        int retVal;
        
        //Get the time in UTC to look up the parameters by going to TAI and
        //then UTC.
        retVal=iauTttai(TT1, TT2, &JulUTC[0], &JulUTC[1]);
        if(retVal!=0) {
            mexErrMsgTxt("An error occurred computing TAI.");
        }
        retVal=iauTaiutc(JulUTC[0], JulUTC[1], &JulUTC[0], &JulUTC[1]);
        switch(retVal){
            case 1:
                mexWarnMsgTxt("Dubious Date entered.");
                break;
            case -1:
                mexErrMsgTxt("Unacceptable date entered");
                break;
            default:
                break;
        }
        
        JulUTCMATLAB[0]=doubleMat2Matlab(&JulUTC[0],1,1);
        JulUTCMATLAB[1]=doubleMat2Matlab(&JulUTC[1],1,1);

        //Get the Earth orientation parameters for the given date.
        mexCallMATLAB(2,retVals,2,JulUTCMATLAB,"getEOP");
        mxDestroyArray(JulUTCMATLAB[0]);
        mxDestroyArray(JulUTCMATLAB[1]);
        
        //%We do not need the polar motion coordinates.
        mxDestroyArray(retVals[0]);
        
        checkRealDoubleArray(retVals[1]);
        if(mxGetM(retVals[1])!=2||mxGetN(retVals[1])!=1) {
            mxDestroyArray(retVals[1]);
            mexErrMsgTxt("Error using the getEOP function.");
            return;
        }
        
        dXdY=mxGetDoubles(retVals[1]);
        dX=dXdY[0];
        dY=dXdY[1];
        
        //Free the returned arrays.
        mxDestroyArray(retVals[1]);
    } else {//Get the celestial pole offsets
        size_t dim1, dim2;
        
        checkRealDoubleArray(prhs[4]);
        dim1 = mxGetM(prhs[4]);
        dim2 = mxGetN(prhs[4]);
        
        if((dim1==2&&dim2==1)||(dim1==1&&dim2==2)) {
            double *dXdY=mxGetDoubles(prhs[4]);
        
            dX=dXdY[0];
            dY=dXdY[1];
        } else {
            mexErrMsgTxt("The celestial pole offsets have the wrong dimensionality.");
            return;
        }
    }
    
    {
    double x, y, s;
        
    //Get the X,Y coordinates of the Celestial Intermediate Pole (CIP) and
    //the Celestial Intermediate Origin (CIO) locator s, using the IAU 2006
    //precession and IAU 2000A nutation models.
    iauXys06a(TT1, TT2, &x, &y, &s);
    
    //Add the CIP offsets.
    x += dX;
    y += dY;
    
    //Get the GCRS-to-CIRS matrix
    iauC2ixys(x, y, s, GCRS2CIRS);
    }
    
    //Allocate space for the return vectors.
    retMat=mxCreateDoubleMatrix(numRow,numVec,mxREAL);
    retData=mxGetDoubles(retMat);
    
    {
        size_t curVec;
        for(curVec=0;curVec<numVec;curVec++) {
            //Multiply the position vector with the rotation matrix.
            iauRxp(GCRS2CIRS, xVec+numRow*curVec, retData+numRow*curVec);
            
            //If a velocity vector was given.
            if(numRow>3) {
                double *velGCRS=xVec+numRow*curVec+3;//Velocity in GCRS
                double *retDataVel=retData+numRow*curVec+3;
                
                //Convert velocity from GCRS to CIRS.
                iauRxp(GCRS2CIRS, velGCRS, retDataVel);
            }
        }
    }
    plhs[0]=retMat;
    
    //If the rotation matrix is desired on the output.
    if(nlhs>1) {
        double *elPtr;
        size_t i,j;
        
        plhs[1]=mxCreateDoubleMatrix(3,3,mxREAL);
        elPtr=mxGetDoubles(plhs[1]);
        
        for (i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                elPtr[i+3*j]=GCRS2CIRS[i][j];
            }
        }
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
