/**GETELLIPSHARMAXES A C++ version of a function to find the Cartesian unit
*                    vectors corresponding to the normalized  gradients of
*                    the components of a point in ellipsoidal harmonic
*                    coordinates. Also, get the magnitudes of the
*                    gradients.
*
*INPUTS: pointHarmon A point given in terms of ellipsoidal harmonic reduced
*                    latitude and longitude in radians and a semi-major
*                    axis in meters.
*                  E The linear eccentricity defining the ellipsoidal
*                    harmonic coordinate system. If this parameter is
*                    omitted, then the linear eccentricity of the WGS84
*                    reference ellipsoid is used.
*
*OUTPUTS: u u(:,1), u(:,2) and u(:,3) are the unit vectors in the
*           respective directions of the respective gradients of reduced
*           latitude, longitude and the semi-major axis.
*         c c(1), c(2) and c(3) are the respective magnitudes of the
*           derivative of the Cartesian position with respect to reduced
*           latitude, longitude and the semi-major axis.
*
*The ellipsoidal harmonic coordinate system is described in Chapter 1.15
*of [1] and is an orthogonal coordinate system.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*[u,c]=getEllipsHarmAxes(pointHarmon,E);
*
*REFERENCES:
*[] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
*SpringerWienNewYork, 2006.
*
*February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "CoordFuncs.hpp"
#include <math.h>
        
void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double E, *pointHarmon, *u,*c;
    mxArray *uMATLAB, *cMATLAB;
    
    if(nrhs<1||nrhs>2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>2) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    if(mxGetM(prhs[0])!=3||mxGetN(prhs[0])!=1) {
       mexErrMsgTxt("The point has the wrong dimensionality.");
       return;
    }
    
    checkRealDoubleArray(prhs[0]);
    pointHarmon=mxGetDoubles(prhs[0]);
    
    //If the linear eccentricity is not given, then use the WGS-84 value
    //from the Constants class.
    if(nrhs<2&&!mxIsEmpty(prhs[1])) {
        mxArray *constantClass, *aMATLAB, *fMATLAB;
        double a,b,f;
        
        mexCallMATLAB(1,&constantClass,0,NULL,"Constants");//Load the Constants class.
        aMATLAB=mxGetProperty(constantClass, 0,"WGS84SemiMajorAxis");

        if(aMATLAB==NULL) {
            mexErrMsgTxt("A necessary default value is missing from the Constants class.");
        }
        fMATLAB=mxGetProperty(constantClass, 0,"WGS84Flattening");
        
        if(fMATLAB==NULL) {
            mexErrMsgTxt("A necessary default value is missing from the Constants class.");
        }
        
        mxDestroyArray(constantClass);
        a=getDoubleFromMatlab(aMATLAB);
        f=getDoubleFromMatlab(aMATLAB);
        b=a*(1 - f);
        E=sqrt(a*a-b*b);
    } else {
        E=getDoubleFromMatlab(prhs[1]);
    }
    
    //Allocate the return variables.
    uMATLAB=mxCreateDoubleMatrix(3, 3,mxREAL);
    cMATLAB=mxCreateDoubleMatrix(3, 1,mxREAL);
    u=mxGetDoubles(uMATLAB);
    c=mxGetDoubles(cMATLAB);
    
    //Compute the values to return.
    getEllipsHarmAxesCPP(u,c,pointHarmon,E);
    
    //Set the return values.
    plhs[0]=uMATLAB;
    if(nlhs>1) {
        plhs[1]=cMATLAB;
    } else {
        mxDestroyArray(cMATLAB);
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

