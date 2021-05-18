/**ITRS2TIRS Rotate vectors from the International terrestrial Reference
*            System (ITRS) into the  Terrestrial Intermediate Reference
*            System (TIRS). The ITRS is essentially the WGS-84 coordinate
*            system: it defines locations with respect to the crust of a
*            non-rotating Earth, where the z axis passes through a fixed
*            point on the surface. On the other hand, the TIRS is nearly
*            the same except the z axis is the axis of rotation of the
*            Earth, which slowly varies over time. Note that the velocity
*            conversion does not include the (small) centrifugal effect of
*            polar motion.
*
*INPUTS: x The NXnumVec collection of vectors in TIRS coordinates to
*          convert (units do not matter). N can be 3, or 6. If the vectors
*          are 3D, then they are position. 6D vectors are assumed to be
*          position and velocity. Since the TIRS and ITRS co-rotate, there
*          is no Coriolis effect to add. Also, the accelerations due to the
*          wobble of the rotation axis over time are not considered. These
*          accelerations are very small. Thus, the function just rotates
*          both halves of the vector.
* TT1, TT2 Two parts of a Julian date given in terrestrial time (TT).
*          The units of the date are days. The full date is the sum of
*          both terms. The date is broken into two parts to provide more
*          bits of precision. It does not matter how the date is split.
*     xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
*          including the effects of tides and librations. If this
*          parameter is omitted or if an empty matrix is passed, the value
*          from the function getEOP will be used.
*
*OUTPUTS: vITRS The NXnumVec vector of values of x rotated from the ITRS
*               into the TIRS.
*        rotMat The 3X3 rotation matrix used to rotate vectors from the
*               ITRS into the TIRS.
*
*The conversion functions from the International Astronomical Union's
*(IAU) Standard's of Fundamental Astronomy library are put together to get
*the necessary rotation matrix.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*[vITRS,rotMat]=ITRS2TIRS(vTIRS,TT1,TT2)
*or if more parameters are known, using the format
*[vITRS,rotMat]=ITRS2TIRS(vTIRS,TT1,TT2,xpyp);
*
*Different celestial coordinate systems are compared in [1].
*
*REFERENCES:
*[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
*    Temporal Coordinate Systems for Target Tracking," Formal Report,
*    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
*    173 pages.
*
*March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t numRow,numVec;
    mxArray *retMat;
    double *retData;
    double ITRS2TIRS[3][3];//Inverse polar motion matrix
    double *xVec, TT1, TT2;
    double xp=0;
    double yp=0;//The polar motion coordinates
    
    if(nrhs<3||nrhs>4){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    numRow=mxGetM(prhs[0]);
    numVec=mxGetN(prhs[0]);
    
    if(!(numRow==3||numRow==6)) {
        mexErrMsgTxt("The input vector has a bad dimensionality.");
    }
    
    checkRealDoubleArray(prhs[0]);
    xVec=mxGetDoubles(prhs[0]);
    
    TT1=getDoubleFromMatlab(prhs[1]);
    TT2=getDoubleFromMatlab(prhs[2]);
    //If xpyp should be found using the function getEOP.
   if(nrhs<4||mxIsEmpty(prhs[3])) {
        mxArray *retVals[1];
        double *xpyp;
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
        mexCallMATLAB(1,retVals,2,JulUTCMATLAB,"getEOP");
        mxDestroyArray(JulUTCMATLAB[0]);
        mxDestroyArray(JulUTCMATLAB[1]);
        
        checkRealDoubleArray(retVals[0]);
        if(mxGetM(retVals[0])!=2||mxGetN(retVals[0])!=1) {
            mxDestroyArray(retVals[0]);
            mexErrMsgTxt("Error using the getEOP function.");
            return;
        }
        
        xpyp=mxGetDoubles(retVals[0]);
        xp=xpyp[0];
        yp=xpyp[1];

        //Free the returned array.
        mxDestroyArray(retVals[0]);
    }
    
     //Get polar motion coordinates, if given.
    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        size_t dim1, dim2;
        
        checkRealDoubleArray(prhs[3]);
        dim1=mxGetM(prhs[3]);
        dim2=mxGetN(prhs[3]);
        
        if((dim1==2&&dim2==1)||(dim1==1&&dim2==2)) {
            double *xpyp=mxGetDoubles(prhs[3]);
        
            xp=xpyp[0];
            yp=xpyp[1];
        } else {
            mexErrMsgTxt("The polar motion coordinates have the wrong dimensionality.");
            return;
        }
    }

    //Get the rotation matrix from TIRS to ITRS.
    {
        double sp;
        double TIRS2ITRS[3][3];//Polar motion matrix
        //Get the Terrestrial Intermediate Origin (TIO) locator s' in
        //radians
        sp=iauSp00(TT1,TT2);
        
        //Get the polar motion matrix
        iauPom00(xp,yp,sp,TIRS2ITRS);
        
        //The inverse polar motion matrix is given by the transpose of the
        //polar motion matrix.
        iauTr(TIRS2ITRS, ITRS2TIRS); 
    }
    
    //Allocate space for the return vectors.
    retMat=mxCreateDoubleMatrix(numRow,numVec,mxREAL);
    retData=mxGetDoubles(retMat);

    {
    size_t curVec;
    for(curVec=0;curVec<numVec;curVec++) {
        //Multiply the position vector with the rotation matrix.
        iauRxp(ITRS2TIRS, xVec+numRow*curVec, retData+numRow*curVec);
        
        //Multiply the velocity vector with the rotation matrix.
        if(numRow>3) {
            double *velITRS=xVec+numRow*curVec+3;//Velocity in ITRS
            double *retDataVel=retData+numRow*curVec+3;//Velocity in ITRS

            //Convert velocity from ITRS to TIRS.
            iauRxp(ITRS2TIRS, velITRS, retDataVel);
        }
    }
    }

    plhs[0]=retMat;
    if(nlhs>1) {
        double *elPtr;
        size_t i,j;
        
        plhs[1]=mxCreateDoubleMatrix(3,3,mxREAL);
        elPtr=mxGetDoubles(plhs[1]);
        
        for (i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                elPtr[i+3*j]=ITRS2TIRS[i][j];
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
