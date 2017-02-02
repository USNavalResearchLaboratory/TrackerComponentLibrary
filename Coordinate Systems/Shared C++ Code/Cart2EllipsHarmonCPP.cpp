/*CART2ELLIPSHARMONCPP  A C++ function to convert from Cartesian to
 *                      ellipsoidal harmonic coefficients.
 *
 *INPUTS: pointsHarmon A pointer to an array of doubles of length
 *                     3*numPoints to hold the converted points with
 *                     elements in the order of
 *                     [reduced latitude;longitude;semiminor axis]
 *        cartPoints   A pointer to an array of 3*numPoints doubles holding
 *                     [x;y;z] values of each of the Cartesian points to
 *                     convert.
 *        numPoints    The number of points that are to be converted.
 *                E    The linear eccentricity of the reference ellipsoid
 *                    that defines the ellipsoid harmonic system.
 *
 *OUTPUTS: None. The results are placed in pointsHarmon.
 *
 *Further comments are given in the Matlab file cart2EllipsHarmon.
 *
 *February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
#include <math.h>
#include <stddef.h>
#include <limits>

using namespace std;
void Cart2EllipsHarmonCPP(double *pointsHarmon,const double *cartPoints, const size_t numPoints, const double E) {
    const double pi=3.14159265358979323846264338327950288419716939937510;
    const double E2=E*E;
    const double *curVec;
    double *curRetVec;
    size_t i;

    curVec=cartPoints;
    curRetVec=pointsHarmon;
    for(i=0;i<numPoints;i++) {
        double x,y,z,t,uPossible,phi;
        double *beta=curRetVec;
        double *lambda=curRetVec+1;
        double *u=curRetVec+2;

        x=curVec[0];
        y=curVec[1];
        z=curVec[2];

        t=x*x+y*y+z*z-E2;
        {
            double innerTerm=(1.0+sqrt(1.0+4.0*E2*z*z/(t*t)));
            //This deals with the case where t is very close to zero so
            //division by t can be problematic.
            if(innerTerm==numeric_limits<double>::infinity()) {
                uPossible=0;
            } else {
                uPossible=sqrt((fabs(t)/2.0)*innerTerm);
            }
        }
        
        *lambda=atan2(y,x);

        if(t>0) {
            double rat;
            *u=uPossible;
            rat=z/(*u);
            
            //Deal with finite precision issues near the poles.
            if(fabs(rat)>1.0) {
                if(z>=0){
                    rat=1.0;
                }else{
                    rat=-1.0;
                }
            }
            
            phi=acos(rat);
        }else {
            double signVal;
            
            if(z>0){
                signVal=1.0;
            }else{
                signVal=-1.0;
            }
            
            phi=acos(signVal*uPossible/E);
            *u=z/cos(phi);
        }

        *beta=pi/2.0-phi;

        curVec+=3;
        curRetVec+=3;
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
