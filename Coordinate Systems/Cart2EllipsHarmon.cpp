/**CART2ELLIPSHARMON A C++ version of a function to convert from Cartesian
*                    coordinates to ellipsoidal harmonic coordinates. Note
*                    that ellipsoidal harmonic coordinates are not the same
*                    as ellipsoidal coordinates.
*
*INPUTS: cartPoints A matrix of the points in ECEF Cartesian coordinates
*                   that are to be transformed into ellipsoidal
*                   coordinates. Each columns of cartPoints is of the
*                   format [x;y;z].
*                 E The linear eccentricity defining the ellipsoidal
*                   harmonic coordinate system. If this parameter is
*                   omitted, then the linear eccentricity of the WGS84
*                   reference ellipsoid is used.
*
*OUTPUTS: pointsHarmon A matrix of the converted points. Each column of the
*                   matrix has the format
 *                  [reduced latitude;longitude;semiminor axis], with
*                   reduced latitude and longitude given in radians.
*
*The ellipsoidal harmonic coordinate system is described in Chapter 1.15 of
*[1]. Note that some folks use the complement of the reduced latitude in
*place of the reduced latitude. The complement is pi/2-the reduced
*latitude.
*
*The conversion should work for all points that are not at the origin. For
*points that are particularly close to the origin, (generally, points deep
*within the Earth) the term
*t=x.^2+y.^2+z.^2-E^2;
*will be negative. If t>0, it can be easily verified that
*u=sqrt((t/2)*(1+sqrt(1+4*E^2*z^2/t^2)));
*On the other hand, for t<0, it can be verified (using, for example,
*Mathematica) that the term
*sqrt((abs(t)/2).*(1+sqrt(1+4*E^2*z.^2./t.^2)));
*is equal to E*abs(cos(phi)). Thus, a different conversion is used
*depending on the value of t. The conversion follows from solving the
*equations in the aforementioned book by Hofmann-Wellenhof and Moritz.
*
*The algorithm can be compiled for use in Matlab using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*pointsHarmon=Cart2EllipsHarmon(cartPoints,E);
*
*REFERENCES:
*[1] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
*    SpringerWienNewYork, 2006.
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
    size_t numPoints;
    mxArray *pointsHarmonMATLAB;
    double E, *pointsHarmon, *cartPoints;
    
    if(nrhs<1||nrhs>2) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    numPoints=mxGetN(prhs[0]);
    if(mxGetM(prhs[0])!=3||numPoints<1) {
       mexErrMsgTxt("The point has the wrong dimensionality.");
       return;
    }

    checkRealDoubleArray(prhs[0]);
    cartPoints=mxGetDoubles(prhs[0]);
    
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
    pointsHarmonMATLAB=mxCreateDoubleMatrix(3, numPoints,mxREAL);
    pointsHarmon=mxGetDoubles(pointsHarmonMATLAB);
    
    //Compute the values to return.
    Cart2EllipsHarmonCPP(pointsHarmon,cartPoints, numPoints, E);

    //Set the return values.
    plhs[0]=pointsHarmonMATLAB;
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
