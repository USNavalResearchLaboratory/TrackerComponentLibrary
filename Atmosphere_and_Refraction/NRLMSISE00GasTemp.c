/**NRLMSISE00GASTEMP  Obtain the gas constituents and temperature
 *                    parameters of the atmosphere at a particular location
 *                    and time using the NRLMSISE-00 atmospheric model for
 *                    dry air.
 *
 *INPUTS: dayOfYear The integer day of the year in the Gregorian calendar
 *                  in universal coordinated time (UTC). Counting starts at
 *                  1. the resolution of the model is not sufficient for it
 *                  to matter whether 365 or 366 is given at the day if it
 *                  isn't/is a leap year.
 *   secondOfTheDay The second of the day. This starts at zero. The
 *                  resolution of the model is not high enough for leap
 *                  seconds to matter, so values above 86400.0 are just
 *                  clipped to 86400.0.
 *        latLonAlt The location under consideration given in WGS-84
 *                  ellipsoidal coordinates of latitude and longitude in
 *                  radians and ellipsoidal height in meters.
 *               Ap An optional parameter specifying the average daily
 *                  geomagnetic index. If omitted, the value 4.0, which is
 *                  suitable for models below 90km is used. The index at a
 *                  particular time can be obtained from a service of the
 *                  International Service of Geomagnetic Indices at
 *                  http://www-app3.gfz-potsdam.de/kp_index/
 *                  A scalar value of AP indicates standard mode. If a 7X1
 *                  vector is passed, then it is assumed that  storm mode
 *                  is desired. In storm mode, Ap(1) is the daily Ap index.
 *                  Ap(2) through Ap(5) are 3 hours Ap indices respectively
 *                  for the current time and 3, 6, and 9 hours before the
 *                  current time. Ap(6) is the average of 8 3-hour Ap
 *                  indices from 12 to 33 hours prior to the current time
 *                  and Ap(7) is the average of 8 3-hour Ap indices from 36
 *                  to 59 hours prior to the current time.
 *             F107 An optional parameter specifying the 10.7cm solar radio
 *                  flux for the previous day. If omitted, a value of 150
 *                  is used, which should be suitable for models below
 *                  90km. Values of the index can be obtained from
 *                  http://www.swpc.noaa.gov/ftpdir/indices/quar_DSD.txt
 *            F107A An optional parameter specifying the 81 day average of
 *                  the 10.7cm radio solar flux. If omitted, a value of 150
 *                  is used, which should be suitable for models below
 *                  90km.
 *
 *OUTPUTS: gasTable A cell array of constituent elements of the atmosphere
 *                  and their number densities. gasTable{i,1} is a string
 *                  listing the name of the ith constituent element of the
 *                  atmosphere. This can take the values
 *                  'He' Helium
 *                  'O'  Elemental Oxygen
 *                  'N2' Nitrogen
 *                  'O2' Oxygen
 *                  'Ar' Argon
 *                  'H'  Elemental Hydrogen
 *                  'N'  Elemental Nitrogen
 *                  'O*' Anomalous Oxygen
 *                  gasTable{i,2} is the corresponding number density of
 *                  the ith constituent element in units of atoms per cubic
 *                  meter.
 *                t Temperature variables. t(1) is the exospheric
 *                  temperature and t(2) is the temperature at altitude.
 *                  The temperature is in units of Kelvin.
 *                d A 9X1 vector containing the same information as
 *                  gasTable but without labels. The entries in d are
 *                  d(1) - HE NUMBER DENSITY
 *                  d(2) - O NUMBER DENSITY
 *                  d(3) - N2 NUMBER DENSITY
 *                  d(4) - O2 NUMBER DENSITY
 *                  d(5) - AR NUMBER DENSITY                       
 *                  d(6) - TOTAL MASS DENSITY  without d(9)
 *                  d(7) - H NUMBER DENSITY
 *                  d(8) - N NUMBER DENSITY
 *                  d(9) - Anomalous oxygen NUMBER DENSITY(M-3)
 *                  where all number densities are in units of atoms per
 *                  cubic meter and the units of d(6) are kilograms per
 *                  cubic meter. Note that d(6) is computed using a
 *                  definition of the atomic mass unit and definitions of
 *                  the atomic masses that are not necessarily the most up
 *                  to date. 
 *
 *This function is essentially a wrapper for the C-code implementation of
 *the  NRLMSISE-00 atmospheric model using appropriate functions to put the
 *time and location in the correct format. The model is described in [1].
 *
 *The algorithm can be compiled for use in Matlab  using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt);
 *or
 *[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP);
 *or
 *[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP,F107);
 *or
 *[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP,F107,F107A);
 *
 *EXAMPLE 1:
 *dayOfYear=172;
 *secondOfTheDay=29000;
 *latLonAlt=[60*(pi/180);-70*(pi/180);400e3];
 *F107=150;
 *F107A=150;
 *AP=4;
 *[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP,F107,F107A);
 *
 *EXAMPLE 2
 *In most cases, one might be given a Julian Date or a Calendar Date.
 *Here, we convert the date appropriately.
 *year=2012;
 *month=11;
 *day=16;
 *hour=0;
 *minute=0;
 *second=0;
 *F107=150;
 *F107A=150;
 *AP=4;
 *latLonAlt=[60*(pi/180);-70*(pi/180);400e3];
 *[Jul1,Jul2]=Cal2UTC(year,month,day,hour,minute,second);
 *[~,dayOfYear,secondOfTheDay]=UTC2DayCount(Jul1,Jul2);
 *[gasTable,t,d]=NRLMSISE00GasTemp(dayOfYear,secondOfTheDay,latLonAlt,AP,F107,F107A);
 *
 *REFERENCES:
 *[1] J. M. Picone, A. E. Hedin, D. P. Drob, and A. C. Aikin, "NRLMSISE-00
 *    empirical model of the atmosphere: Statistical comparisons and
 *    scientific issues," Journal of Geophysical Research: Space Physics,
 *    vol. 107, no. A12, Dec. 2002.
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

 /*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"
#include "nrlmsise-00.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //Multiplicative constant for converting from radians to degrees.
    const double rad2Deg=57.295779513082320876798154814105;
    struct nrlmsise_flags theFlags;
    struct nrlmsise_output theOutput;
    struct nrlmsise_input theInput;
    struct ap_array theAPData;
    double Jul1,Jul2, *latLonAlt,ap,f107,f107A;
    size_t i;

    //This activates everything and makes the output in meters.
    for(i=0;i<24;i++) {
        theFlags.switches[i]=1;
    }
    
    if(nrhs<3||nrhs>6){
        mexErrMsgTxt("Wrong number of inputs.");
    }
    
    if(nlhs>3){
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    theInput.doy=getIntFromMatlab(prhs[0]);
    theInput.sec=getDoubleFromMatlab(prhs[1]);
    if(theInput.sec>86400.0) {
        theInput.sec=86400.0;
    }
    //The year argument is not used by the gtd7 function.
    theInput.year=0;
    
    checkRealDoubleArray(prhs[2]);
    if(mxGetM(prhs[2])!=3||mxGetN(prhs[2])!=1) {
        mexErrMsgTxt("The latitude, longitude, and altitude point has the wrong dimensionality.");
    }
    latLonAlt=mxGetDoubles(prhs[2]);
    
    ap=4;//The default value for the geomagnetic activity index.
    
    //The default values for the 10.7 cm radio flux parameters.
    f107=150;
    f107A=150;
    
    if(nlhs>3&&!mxIsEmpty(prhs[3])) {
        checkRealDoubleArray(prhs[3]);
        if(mxGetM(prhs[3])==1&&mxGetM(prhs[3])==1) {
            ap=getDoubleFromMatlab(prhs[3]);
        } else if(mxGetM(prhs[3])==7&&mxGetM(prhs[3])==1) {
            double *aVec=mxGetDoubles(prhs[3]);
            for(i=0;i<7;i++) {
                theAPData.a[i]=aVec[i];
            }

            ap=theAPData.a[0];
            theFlags.switches[9]=-1;//This activates storm mode, indicating the presence of the ap_array structure. 
        } else {
            mexErrMsgTxt("The AP parameter has the wrong dimensionality.");
        }
    }
    
    if(nlhs>4&&!mxIsEmpty(prhs[4])) {
        f107=getDoubleFromMatlab(prhs[4]);
    }
    
    if(nlhs>5&&!mxIsEmpty(prhs[5])) {
        f107=getDoubleFromMatlab(prhs[5]);
    }
    
    theInput.f107A=f107A;
    theInput.f107=f107;
    theInput.ap=ap;
    theInput.ap_a=&theAPData;
    
    //The 1000 transforms the units into kilometers as required by the
    //gtd7 function.
    theInput.alt=latLonAlt[2]/1000.0;

    //The 180/pi transforms the units from radians to degrees, as
    //required by the gtd7 function.
    theInput.g_lat=latLonAlt[0]*rad2Deg;
    theInput.g_long=latLonAlt[1]*rad2Deg;
    
    //This is deterministically set as described in THE NRLMSISE-00 AND
    //HWM-93 USERS GUIDE: Version 1.50, November 2003 by Douglas P. Drob.
    //Basically, the precision of the model is not high enough for it to
    //matter if the Equation of Time is used as opposed to a simple
    //approximation of local solar time.
    theInput.lst=16;//theInput.sec/3600+theInput.g_long/15;

    ///Run the routine
    gtd7(&theInput,&theFlags, &theOutput);
    
    //Extract the output to return.
    {
        mxArray *celArray2Ret=mxCreateCellMatrix(8,2);
        
        //Put all of the names and mass densitites of the constituent
        //elements into the cell array.
        mxSetCell(celArray2Ret,0,mxCreateString("He"));
        mxSetCell(celArray2Ret,8,doubleMat2Matlab(&theOutput.d[0],1,1));
        mxSetCell(celArray2Ret,1,mxCreateString("O"));
        mxSetCell(celArray2Ret,9,doubleMat2Matlab(&theOutput.d[1],1,1));
        mxSetCell(celArray2Ret,2,mxCreateString("N2"));
        mxSetCell(celArray2Ret,10,doubleMat2Matlab(&theOutput.d[2],1,1));
        mxSetCell(celArray2Ret,3,mxCreateString("O2"));
        mxSetCell(celArray2Ret,11,doubleMat2Matlab(&theOutput.d[3],1,1));
        mxSetCell(celArray2Ret,4,mxCreateString("Ar"));
        mxSetCell(celArray2Ret,12,doubleMat2Matlab(&theOutput.d[4],1,1));
        mxSetCell(celArray2Ret,5,mxCreateString("H"));
        mxSetCell(celArray2Ret,13,doubleMat2Matlab(&theOutput.d[6],1,1));
        mxSetCell(celArray2Ret,6,mxCreateString("N"));
        mxSetCell(celArray2Ret,14,doubleMat2Matlab(&theOutput.d[7],1,1));
        mxSetCell(celArray2Ret,7,mxCreateString("O*"));
        mxSetCell(celArray2Ret,15,doubleMat2Matlab(&theOutput.d[8],1,1));
        
        plhs[0]=celArray2Ret;
        if(nlhs>1) {
            plhs[1]=doubleMat2Matlab(theOutput.t,2, 1);
            
            if(nlhs>2) {
                plhs[2]=doubleMat2Matlab(theOutput.d,9, 1);
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
