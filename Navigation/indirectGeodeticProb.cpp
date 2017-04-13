/**INDIRECTGEODETICPROB  Solve the indirect geodetic problem. That is, given
 *                       two points on an ellipsoidal Earth, find the
 *                       initial bearing and distance one must travel to
 *                       take the shortest (geodesic) path between the
 *                       points.
 *
 *INPUTS:latLonStart A 2XN matrix of the N initial points given in geodetic
 *                  latitude and longitude in radians of the format
 *                  [latitude;longitude] for each column (point). The
 *                  latitude must be between -pi/2 and pi/2 radians and the
 *                  longitude between -pi and pi radians.
 *        latLonEnd The 2XN matrix of final points for the geodesic path
 *                  given in geodetic latitude  and longitude in radians.
 *                  latLonEnd has the same format at latLonStart.
 *                a The semi-major axis of the reference ellipsoid (in
 *                  meters). If this argument is omitted, the value in
 *                  Constants.WGS84SemiMajorAxis is used.
 *                f The flattening factor of the reference ellipsoid. If
 *                  this argument is omitted, the value in
 *                  Constants.WGS84Flattening is used.
 *
 *OUTPUTS: azStart  The NX1 forward azimuth at the starting points in 
 *                  radians East of true North on the reference ellipsoid.
 *                  This is the initial heading one would travel to go
 *                  between latLonStart and latLonEnd for each of the point
 *                  pairs
 *             dist The NX1 vector of geodetic distances between the
 *                  starting and stopping points in meters.
 *            azEnd The forward azimuth at the ending point in radians
 *                  East of true North on the reference ellipsoid.
 *
 *The function is essentially a Matlab interface for the implementation in
 *GeographicLib, which is documented in [1], [2], and [3]. GeographicLib
 *can be downloaded from
 *http://geographiclib.sourceforge.net
 *Though a native Matlab version of the relevant function in GepgraphicLib
 *exists, it is rather slow. hence the need for this interface to the
 *compiled version.
 *
 *The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[azStart,dist,azEnd]=indirectGeodeticProb(latLonStart,latLonEnd);
 *or if something other than the WGS84 reference ellipsoid is used
 *[azStart,dist,azEnd]=indirectGeodeticProb(latLonStart,latLonEnd,a,f);
 *
 *REFERENCES:
 *[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
 *    87, no. 1, pp. 43?45, Jan. 2013. [Online].
 *    Available: http://arxiv.org/pdf/1109.4448.pdf
 *[2] C. F. F. Karney. (2013, 31 Aug.) Addenda and errata for papers on
 *    geodesics. [Online].
 *    Available: http://geographiclib.sourceforge.net/geod-addenda.html
 *[3] C. F. F. Karney. (2011, 7 Feb.) Geodesics on an ellipsoid of
 *    revolution.
 *    [Online]. Available: http://arxiv.org/pdf/1102.1215.pdf
 *
 *April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
 */
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
#include "MexValidation.h"

//For the wrapRange function so that there is no issue if latitudes and
//longitudes are given outside of the valid range for the GeographicLib
//functions.
#include "mathFuncs.hpp"

//For fabs
#include <cmath>

#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicExact.hpp>

using namespace GeographicLib;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const double pi=3.1415926535897932384626433832795;
    size_t numPoints;
    double *latLonStart, *latLonEnd;
    double a,f;
    double *azStart, *dist, *azEnd;
    mxArray *azStartMATLAB, *distMATLAB, *azEndMATLAB;
    
    if(nrhs<2||nrhs>4){
        mexErrMsgTxt("Wrong number of inputs");
    }

    if(nlhs>3) {
        mexErrMsgTxt("Wrong number of outputs.");
    }
    
    checkRealDoubleArray(prhs[0]);
    checkRealDoubleArray(prhs[1]);
    
    {
        size_t numRow;
        numRow = mxGetM(prhs[0]);
        numPoints = mxGetN(prhs[0]);
        
        if(numRow!=2) {
            mexErrMsgTxt("The latLonStart vector has a bad dimensionality.");
        }
        
        numRow = mxGetM(prhs[1]);
        
        if(numRow!=2||mxGetN(prhs[1])!=numPoints) {
            mexErrMsgTxt("The latLonEnd vector has a bad dimensionality.");
        }
    }
    
    //Get the starting and ending points.
    latLonStart=reinterpret_cast<double*>(mxGetData(prhs[0]));
    latLonEnd=reinterpret_cast<double*>(mxGetData(prhs[1]));
    
    if(nrhs>2) {
        a=getDoubleFromMatlab(prhs[2]);
    } else {
        a=getScalarMatlabClassConst("Constants","WGS84SemiMajorAxis");
    }

    if(nrhs>3) {
        f=getDoubleFromMatlab(prhs[3]);
    } else {
        f=getScalarMatlabClassConst("Constants","WGS84Flattening");
    }
    
    //Allocate space for the return variables.
    azStartMATLAB=mxCreateDoubleMatrix(numPoints,1,mxREAL);
    distMATLAB=mxCreateDoubleMatrix(numPoints,1,mxREAL);
    azEndMATLAB=mxCreateDoubleMatrix(numPoints,1,mxREAL);
    
    azStart=reinterpret_cast<double*>(mxGetData(azStartMATLAB));
    dist=reinterpret_cast<double*>(mxGetData(distMATLAB));
    azEnd=reinterpret_cast<double*>(mxGetData(azEndMATLAB));
    
    //Solve the indirect geodetic problem for each of the point pairs.
    //The class used to solve the problem is chosen based on the value of
    //f. As per the documentation is GeographicLib, the Geodesic class is
    //best suited for problems where |f|<=0.01. For other problems, the
    //GeodesicExact class is the better choice.
    
    if(fabs(f)<=0.01) {
        Geodesic geod(a, f);
        size_t curPoint,curBaseIdx;
        
        curBaseIdx=0;
        for(curPoint=0;curPoint<numPoints;curPoint++) {
            //The function takes the values in degrees.
            const double lat1=wrapRangeMirrorCPP(latLonStart[curBaseIdx]*(180.0/pi),-90,90);
            const double lon1=wrapRangeCPP(latLonStart[curBaseIdx+1]*(180.0/pi),-180,180);
            const double lat2=wrapRangeMirrorCPP(latLonEnd[curBaseIdx]*(180.0/pi),-90,90);
            const double lon2=wrapRangeCPP(latLonEnd[curBaseIdx+1]*(180.0/pi),-180,180);
            double s12,azi1,azi2;//Return values are put in here.
            
            geod.Inverse(lat1, lon1, lat2, lon2, s12, azi1, azi2);
            
            azStart[curPoint]=azi1*(pi/180.0);
            dist[curPoint]=s12;
            azEnd[curPoint]=azi2*(pi/180.0);
            
            curBaseIdx+=2;
        }
    } else {//For the case where |f|>0.01. All that changes is the class
            //used.
        GeodesicExact geod(a, f);
        size_t curPoint,curBaseIdx;
        
        curBaseIdx=0;
        for(curPoint=0;curPoint<numPoints;curPoint++) {
            //The function takes the values in degrees.
            const double lat1=latLonStart[curBaseIdx]*(180.0/pi);
            const double lon1=latLonStart[curBaseIdx+1]*(180.0/pi);
            const double lat2=latLonEnd[curBaseIdx]*(180.0/pi);
            const double lon2=latLonEnd[curBaseIdx+1]*(180.0/pi);
            double s12,azi1,azi2;//Return values are put in here.
            
            geod.Inverse(lat1, lon1, lat2, lon2, s12, azi1, azi2);
            
            azStart[curPoint]=azi1*(pi/180.0);
            dist[curPoint]=s12;
            azEnd[curPoint]=azi2*(pi/180.0);
            
            curBaseIdx+=2;
        }
    }
    
    //Set the return values.
    plhs[0]=azStartMATLAB;
    
    if(nlhs>1) {
        plhs[1]=distMATLAB;
        
        if(nlhs>2) {
            plhs[2]=azEndMATLAB;
        } else {
            mxDestroyArray(azEndMATLAB);  
        }
    } else {
        mxDestroyArray(distMATLAB);
        mxDestroyArray(azEndMATLAB);
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
