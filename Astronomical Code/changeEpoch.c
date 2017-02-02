/**CHANGEEPOCH Given star catalog data for a particular epoch, change it to
 *             a new epoch. This function is useful for transforming data
 *             in star catalogs to the J2000.0 epoch so that the function
 *             starCat2Obs can be used to determine what an observer on the
 *             Earth should see.
 *
 *INPUTS: catData   Cat data is a matrix of stars (one per row) that are to
 *                  be converted, where each row has the following format:
 *catData(:,1) RArad    Right Ascension in (rad) ICRS at the J2000.0 epoch
 *catData(:,2) DErad    Declination in (rad) ICRS at the J2000.0 epoch
 *catData(:,3) Plx      Parallax (rad)
 *catData(:,4) pmRA     Proper motion in Right Ascension (rad/yr)
 *                      in the form dRA/dt and not cos(Dec)*dRA/dt.
 *catData(:,5) pmDE     Proper motion in Declination (rad/yr)
 *catData(:,6) vRad     Radial velocity of the star in meters per second
 *                      with a positive sign meaning that a star is
 *                      receding (going away from) Earth.
 *Jul1,Jul2 Two parts of a Julian date given in TDB specifying the date of
 *          the original epoch. The units of the date are days. The full
 *          date is the sum of both terms. The date is broken into two
 *          parts to provide more bits of precision. It does not matter how
 *          the date is split. If these parameters are omitted or if empty
 *          matrices are passed, then a source epoch of Jul1=2448349;
 *          Jul2=0.062500019022333; is used. These values are what one
 *          obtains by putting the terrestrial time 2448349.0625, the epoch
 *          of the Hipparcos catalog, though the function TT2TDB.
 *Jul3,Jul4 Two parts of a Julian date for the epoch into which the
 *          parameters should be converted in TDB. If these parameters are
 *          omitted, then a terrestrial time of Jul1=2451545.0; Jul2=0;
 *          will be used, corresponding to the J2000.0 epoch.
 *
 *OUTPUTS:  catData The same as the input catData, except the epoch has
 *                  been changed.
 *
 *This function is pretty much a wrapper for the iauPmsafe that is part of
 *the International Astronomical Union's Standards of Fundamental Astronomy
 *library.
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *catData=changeEpoch(catData,Jul1,Jul2,Jul3,Jul4);
 *
 *March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
/*This header is for the SOFA library.*/
#include "sofa.h"
#include "MexValidation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const double pi=3.1415926535897932384626433832795;
    //The constant to convert arcseconds to radians. The multiplications go
    //as->arcminutes->deg->rad
    const double as2Rad=(1.0/60.0)*(1.0/60.0)*(pi/180.0);
    const double rad2as=1.0/as2Rad;
    double *catData;//The inputs
    mxArray *catDataRetMATLAB;
    double *catDataRet;//The outputs
    //Pointers to the rows of catData
    double *RArad, *DErad, *Plx, *pmRA, *pmDE,*vRad;
    //Pointers to the rows of catDataRet
    double *RArad2, *DErad2, *Plx2, *pmRA2, *pmDE2,*vRad2;
    double Jul1,Jul2, Jul3,Jul4;
    size_t numStars,i;

    if(nrhs<1||nrhs>5||nrhs==2||nrhs==4) {
        mexErrMsgTxt("Wrong number of inputs.");
        return;
    }
    
    if(nlhs>4) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    numStars=mxGetM(prhs[0]);
    if(numStars<1||mxGetN(prhs[0])!=6) {
        mexErrMsgTxt("The catalog data has the wrong dimensionality.");
        return;
    }
    checkRealDoubleArray(prhs[0]);
    catData=(double*)mxGetData(prhs[0]);

    //If non-empty values are provided.
    if(nrhs>1&&mxGetM(prhs[1])!=0) {
        Jul1=getDoubleFromMatlab(prhs[1]);
        Jul2=getDoubleFromMatlab(prhs[2]);
    } else {
         Jul1=2448349.0;
         Jul2=0.062500019022333;
   }
    
    if(nrhs>3) {
        Jul3=getDoubleFromMatlab(prhs[3]);
        Jul4=getDoubleFromMatlab(prhs[4]);
    } else {
         Jul3=2451545.0;
         Jul4=0;
     }
    
    //Allocate the return matrices
    catDataRetMATLAB=mxCreateDoubleMatrix(numStars,6,mxREAL);
    catDataRet=(double*)mxGetData(catDataRetMATLAB);
    
    //Set the pointers for each column of data in the original catalog.
    RArad=catData;
    DErad=catData+numStars;
    Plx=catData+2*numStars;
    pmRA=catData+3*numStars;
    pmDE=catData+4*numStars;
    vRad=catData+5*numStars;
    
    //Set the pointers for each column of data in the converted catalog.
    RArad2=catDataRet;
    DErad2=catDataRet+numStars;
    Plx2=catDataRet+2*numStars;
    pmRA2=catDataRet+3*numStars;
    pmDE2=catDataRet+4*numStars;
    vRad2=catDataRet+5*numStars;

     //Convert the parameters for each of the stars.
     for(i=0;i<numStars;i++) {
         int retVal;
        
        //The function wants the parallax to be in milliarcseconds, but the
        //input is in SI units of radians. Thus, we have to convert the
        //units, run the function and convert the result back. Also, the
        //function wants the radial motion to be in units of kilometers per
        //second, so a similar conversion between meters and kilometers has
        //to happen.
        retVal=iauPmsafe(RArad[i], DErad[i], pmRA[i], pmDE[i],Plx[i]*rad2as, vRad[i]/1000.0, Jul1, Jul2, Jul3, Jul4, &RArad2[i], &DErad2[i], &pmRA2[i], &pmDE2[i],&Plx2[i], &vRad2[i]);
               
        //Convert back to radians.
        Plx2[i]*=as2Rad;
        //Convert back to meters.
        vRad2[i]*=1000.0;
        
        if(retVal<0) {
            mxDestroyArray(catDataRetMATLAB);
            mexErrMsgTxt("An error occurred in the iauPmsafe function");
            return;
        }
     }
    
    //Set the outputs
    plhs[0]=catDataRetMATLAB;
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
