/**LIGHTDEFLECTCORR Rotate vectors to account for light-deflection due to
 *             one or more solar system bodies using an approximate formula
 *             to deal with the light-time problem without resorting to ray
 *             tracing.
 *
 *INPUTS: vObsSource A 3XN set of vectors pointing from the observer to the 
 *              objects being observed without corruption due to light 
 *              deflection by the Sun or a planet These are direction 
 *              vectors in the International Celestial Reference System 
 *              (ICRS). Normally, one would expect vObsSource to be unit
 *              vectors.
 *      posObs The 3X1 Cartesian position vector of the observer in the
 *             Barycentric Celestial Reference System (BCRS) at the time of
 *             observation. The units are meters.
 *      MSolar The numBodX1 vector of the masses of the gravitating bodies
 *             that will deflect the light. The masses are solar masses.
 *             That is, the ratio of the mass of the body to the mass of
 *             the Sun. Thus, the input is unitless.
 *       xBody A 6XnumBodies set of position and velocity state vectors in
 *             the BCRS for the defecting bodies at the time of
 *             observation. The units are meters and meters per second.
 * deflecLimit An optional numBodX1 set of deflection limiter parameters in
 *             radians that are phi^2/2 for each body, where phi is the
 *             angular distance between the source and the thing that is to
 *             deflect the light of the source. As phi goes below the
 *             threshold, the deflection applied goes to zero. If this is
 *             omitted, the default value of 1e-10 for all of the inputs
 *             is used.
 *
 *OUTPUTS: uAberr The 3XN set of N vectors rotated to deal with the
 *            deflection due to the specified astronomical bodies.
 *             
 *This function is primarily a wrapper for the approximate light deflection
 *formula used in the iauLdn in the International Astronomical Union's
 *(IAU's) Standards of Fundamental Astronomy (SOFA) library.
 *
 *Example solar mass values and deflection limit values taken from the
 *IAU's documentation for the iauLdn function are
 *Body       Solar Mass     Deflection Limit
 *Sun        1.0            6e-6
 *Jupiter    0.00095435     3e-9
 *Saturn     0.00028574     3e-10
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *vDefl=lightDeflectCorr(vObsSource,posObs,MSolar,xBody);
 *or
 *vDefl=lightDeflectCorr(vObsSource,posObs,MSolar,xBody,deflecLimit);
 *
 *April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
//For constants in the SOFA library.
#include "sofam.h"
#include <limits.h>
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    size_t numSource, numBodies, baseIdx,curBody,curSource;
    double *vObsSource,posObs[3], *MSolar, *xBody, *deflecLimit=NULL;
    iauLDBODY *bodyParam;
    mxArray *vDeflMATLAB;
    double *vDefl;

    if(nrhs<4||nrhs>5) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }

    //Check the inputs
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    checkRealDoubleArray(prhs[2]);
    
    //Check the inputs
    numSource=mxGetN(prhs[0]);

    if(mxGetM(prhs[0])!=3||numSource==0) {
        mexErrMsgTxt("The input vObsSource has the wrong dimensionality.");
    }

    if(mxGetM(prhs[1])!=3||mxGetN(prhs[1])!=1) {
        mexErrMsgTxt("The input posObs has the wrong dimensionality.");
    }

    numBodies=mxGetM(prhs[2]);
    
    if(numBodies>INT_MAX) {
        mexErrMsgTxt("Too many bodies are specified.");
    }
    

    if(numBodies==0||mxGetN(prhs[2])!=1) {
        mexErrMsgTxt("The input MSolar has the wrong dimensionality.");
    }

    if(mxGetM(prhs[3])!=6||mxGetN(prhs[2])!=numBodies) {
        mexErrMsgTxt("The input xBody has the wrong dimensionality.");
    }

    if(nrhs>3) {
        checkRealDoubleArray(prhs[4]); 

        if(mxGetM(prhs[4])!=numBodies||mxGetN(prhs[4])!=1) {
           mexErrMsgTxt("The input deflecLimit has the wrong dimensionality."); 
        }

    } else {
        deflecLimit=NULL;
    }

    vObsSource=mxGetDoubles(prhs[0]);
            
    //Get the observer position and convert from meters to AU.
    {
        double *temp=mxGetDoubles(prhs[1]);
        posObs[0]=temp[0]/DAU;
        posObs[1]=temp[1]/DAU;
        posObs[2]=temp[2]/DAU;
    }
    
    MSolar=mxGetDoubles(prhs[2]);
    //The units have to be converted to AU and AU/Day
    xBody=mxGetDoubles(prhs[3]);
       
    //Allocate space to hold the parameters of the astronomical bodies in a
    //manner suitable for the iauLdn function.
    bodyParam=(iauLDBODY*)mxMalloc(sizeof(iauLDBODY)*numBodies);
             
    baseIdx=0;
    for(curBody=0;curBody<numBodies;curBody++) {
        size_t i;
        
        bodyParam[curBody].bm=MSolar[curBody];
        if(deflecLimit!=NULL) {
            bodyParam[curBody].dl=deflecLimit[curBody];
        } else {
            //A value suitably small for Saturn.
            bodyParam[curBody].dl=3e-10;
        }
        
        //Position
        for(i=0;i<3;i++) {
            //The division converts from meters to AU.
            bodyParam[curBody].pv[0][i]=xBody[baseIdx+i]/DAU;
        }
        
        //Velocity
        for(i=0;i<3;i++) {
            //Convert from meters per second BCRS to AU/ day.
            bodyParam[curBody].pv[1][i]=xBody[3+baseIdx+i]*(1/DAU)*(1/DAYSEC);
        }
        
        baseIdx+=6;
    }
    
    //Allocate space for the return values.
    vDeflMATLAB=mxCreateDoubleMatrix(3,numSource,mxREAL);
    vDefl=mxGetDoubles(vDeflMATLAB);
    
    baseIdx=0;
    for(curSource=0;curSource<numSource;curSource++) {
        double vecMag, sc[3];//Unit vector to the source
        
        //Get a unit direction vector and magnitude to the current source.
        iauPn(vObsSource+baseIdx, &vecMag, sc);
        
        iauLdn((int)numBodies, bodyParam, posObs, sc,vDefl+baseIdx);
        
        //Deal with possibly non-unit magnitudes on the input.
        vDefl[baseIdx]*=vecMag;
        vDefl[baseIdx+1]*=vecMag;
        vDefl[baseIdx+2]*=vecMag;
        
        baseIdx+=3;
    }
    
    //Free temporary memory; set the return value.
    mxFree(bodyParam);
    plhs[0]=vDeflMATLAB;
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
