/**APPROXSOLARSYSVEC Get the approximate position and velocity vector of
 *                 the Earth or of the planets in coordinates aligned with
 *                 the International Celestial Reference System (ICRS).
 *                 The position for everything is heliocentric, and for the
 *                 Earth an additional barycentric (Like the Barycentric
 *                 Celestial Reference System [BCRS], except TDB is used
 *                 and not TCB) vector is returned. Use the function
 *                 solarBodyVec if higher-precision with respect to
 *                 observers near the Earth is desired.
 *
 *INPUTS: Jul1,Jul2 Two parts of a Julian date given in Barycentric
 *                  dynamical time (TDB). However, terrestrial time (TT)
 *                  can generally be substituted as the two timescales are
 *                  relatively close and the approximation is not very high
 *                  fidelity. The units of the date are days. NX1 or 1XN 
 *                  vectors can be passed if one wishes to find the 
 *                  positions and velocities at multiple times. The full 
 *                  date is the sum of both terms. The date is broken into
 *                  two parts to provide more bits of precision. 
 *                  It does not matter how the date is split. The 
 *                  valid range of dates is from 1900AD-2100AD for the 
 *                  Earth (solarBody=0) and 1000AD- 3000AD for everything
 *                  else.
 *        solarBody A parameter specifying the solar body. Possible values
 *                  are:
 *                  0 (The default if omitted) The Earth.
 *                  1 Mercury
 *                  2 Venus
 *                  3 The Earth-Moon barycenter
 *                  4 Mars
 *                  5 Jupiter
 *                  6 Saturn
 *                  7 Uranus
 *                  8 Neptune
 *
 *OUTPUTS: xHelio  A 6XN state vector of the body where the coordinate
 *                 axes are aligned with the ICRS and the origin is the
 *                 Sun. xHelio(1:3,i) is the position for the ith date and
 *                 xHelio(4:6,i) is the velocity. If an error occurred for
 *                 a particular date, zeros are returned for that date.
 *                 The units are meters and meters per second TDB (which is
 *                 close to TT).
 *         xBary   A 6XN state vector of the body where the coordinate
 *                 axes are aligned with the ICRS and the origin is the
 *                 solar system barycenter. This is essentially the BCRS.
 *                 xBary(1:3,i) is the position for the ith date and
 *                 xBary(4:6,i) is the velocity. If an error occurred for
 *                 a particular date, zeros are returned for that date.
 *                 The units are meters and meters per second TDB (which is
 *                 close to TT). xBary is only returned if solarBody=0.
 *                 Otherwise, an empty matrix is returned.
 *        exitFlag A NX1 set of numbers indicating the exit status.
 *                 Possible values for each date are
 *                 0 = OK
 *                 1 = Date outside of valid range.
 *                 2 = Algorithm failed to converge.
 *
 *This function is just a wrapper for the iauEpv00 function in the
 *International Astronomical Union's (IAU) Standard's of Fundamental
 *Astronomy library.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[xHelio,xBary,exitFlag]=approxSolarSysVec(Jul1,Jul2,solarBody);
 *
 *April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

 /*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *Jul1,*Jul2;
    double AU2Meters;
    double *xHelio, *exitFlag,*xBary=NULL;
    mxArray *xHelioMATLAB, *xBaryMATLAB, *exitFlagMATLAB;
    size_t totalVecs,curVec;
    int solarBody=0;

    if(nrhs<2||nrhs>3){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>3){
        mexErrMsgTxt("Wrong number of outputs");
    }
    
    if(nrhs>2) {
        solarBody=getIntFromMatlab(prhs[2]);
    }
    
    if(solarBody>8||solarBody<0) {
        mexErrMsgTxt("Invalid solar body specified.");
    }
    
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    /*Extract the input arguments.*/
    {
        size_t numRow,numCol;
        numRow=mxGetM(prhs[0]);
        numCol=mxGetN(prhs[0]);
        
        totalVecs=numRow*numCol;
    
        numRow=mxGetM(prhs[1]);
        numCol=mxGetN(prhs[1]);
        
        if(numRow*numCol!=totalVecs){
            mexErrMsgTxt("The input dimensionalitites do not match.");
        }
        
        if(totalVecs<=0) {
            mexErrMsgTxt("Empty matrices were passed instead of dates.");
        }
    }
    
    //Get the dates.
    Jul1=mxGetDoubles(prhs[0]);
    Jul2=mxGetDoubles(prhs[1]);
    
    //Get the astronomical unit constant.
    AU2Meters=getScalarMatlabClassConst("Constants", "AstronomialUnit");
    
    //Allocate space for the return values.
    xHelioMATLAB=mxCreateDoubleMatrix(6,totalVecs,mxREAL);
    xHelio=mxGetDoubles(xHelioMATLAB);
    
    if(solarBody==0) {
        xBaryMATLAB=mxCreateDoubleMatrix(6,totalVecs,mxREAL);
        xBary=mxGetDoubles(xBaryMATLAB);
    } else {
        //The approximations for the planets are not offered up in
        //Barycentric coordinates.
        xBaryMATLAB=mxCreateDoubleMatrix(0,0,mxREAL);
    }
        
    exitFlagMATLAB=mxCreateDoubleMatrix(totalVecs,1,mxREAL);
    exitFlag=mxGetDoubles(exitFlagMATLAB);
    
    //Run the functions.
    if(solarBody==0) {
        size_t curIdx=0;
        for(curVec=0;curVec<totalVecs;curVec++) {
            double pvh[2][3];
            double pvb[2][3];
            
            exitFlag[curVec]=iauEpv00(Jul1[curVec], Jul2[curVec], pvh, pvb);
            
            //Put the results in the return matrices and convert the units.
            //The returned distances are in astronomical units. The
            //returned velocities are in astronomical units per Julian
            //day TDB. Note that there is 86400.0 seconds per Julian day.
            xHelio[curIdx+0]=pvh[0][0]*AU2Meters;
            xHelio[curIdx+1]=pvh[0][1]*AU2Meters;
            xHelio[curIdx+2]=pvh[0][2]*AU2Meters;
            xHelio[curIdx+3]=pvh[1][0]*AU2Meters*(1/86400.0);
            xHelio[curIdx+4]=pvh[1][1]*AU2Meters*(1/86400.0);
            xHelio[curIdx+5]=pvh[1][2]*AU2Meters*(1/86400.0);
            
            xBary[curIdx+0]=pvb[0][0]*AU2Meters;
            xBary[curIdx+1]=pvb[0][1]*AU2Meters;
            xBary[curIdx+2]=pvb[0][2]*AU2Meters;
            xBary[curIdx+3]=pvb[1][0]*AU2Meters*(1/86400.0);
            xBary[curIdx+4]=pvb[1][1]*AU2Meters*(1/86400.0);
            xBary[curIdx+5]=pvb[1][2]*AU2Meters*(1/86400.0);
            
            curIdx=curIdx+6;
        }
    } else {
        size_t curIdx=0;
        for(curVec=0;curVec<totalVecs;curVec++) {
            double pvh[2][3];
            
            exitFlag[curVec]=iauPlan94(Jul1[curVec], Jul2[curVec], solarBody, pvh);
            
            //Put the results in the return matrices and convert the units.
            //The returned distances are in astronomical units. The
            //returned velocities are in astronomical units per Julian
            //day TDB. Note that there is 86400.0 seconds per Julian day.
            xHelio[curIdx+0]=pvh[0][0]*AU2Meters;
            xHelio[curIdx+1]=pvh[0][1]*AU2Meters;
            xHelio[curIdx+2]=pvh[0][2]*AU2Meters;
            xHelio[curIdx+3]=pvh[1][0]*AU2Meters*(1/86400.0);
            xHelio[curIdx+4]=pvh[1][1]*AU2Meters*(1/86400.0);
            xHelio[curIdx+5]=pvh[1][2]*AU2Meters*(1/86400.0);

            curIdx=curIdx+6;
        }
    }
    
    //Set the return values.
    plhs[0]=xHelioMATLAB;
    if(nlhs>1) {
        plhs[1]=xBaryMATLAB;
        if(nlhs>2) {
            plhs[2]=exitFlagMATLAB;
        }
        else {
            mxDestroyArray(exitFlagMATLAB);
        }
    }
    else {
        mxDestroyArray(xBaryMATLAB);
        mxDestroyArray(exitFlagMATLAB);
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
