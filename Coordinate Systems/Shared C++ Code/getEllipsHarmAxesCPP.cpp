/*GETELLIPSHARMAXESCPP  A C++ function to compute the basis vectors of an
 *                      ellipsoidal harmonic coordinate system.
 *
 *INPUTS: uRet, c   Pointers to arrays in which the unit vectors for the
 *                  ENU axes and the magnitudes of the unnormalized vectors
 *                  are placed.
 *     pointHarmon  A pointer to an array of doubles having three elements
 *                  ordered [reduced latitude;longitude;semiminor axis].
 *                E    The linear eccentricity of the reference ellipsoid
 *                    that defines the ellipsoid harmonic system.
 *
 *OUTPUTS: None. The results are placed in uRet and c.
 *
 *Further comments are given in the Matlab file getEllipsHarmAxes.
 *
 *February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "CoordFuncs.hpp"
#include <math.h>
#include <stddef.h>

void getEllipsHarmAxesCPP(double *uRet, double *c,const double *pointHarmon,const double E) {
    const double beta=pointHarmon[0];
    const double lambda=pointHarmon[1];
    const double u=pointHarmon[2];
    const double cosBeta=cos(beta);
    const double sinBeta=sin(beta);
    const double cosLambda=cos(lambda);
    const double sinLambda=sin(lambda);
    const double rVal=sqrt(E*E+u*u);
    double *dBeta,*dLambda,*du;//Pointer to the normalized gradient vectors.
    double *c1, *c2, *c3;//Pointers to the magnitudes of the unnormalized gradient vectors.
    size_t i;
    
    dBeta=uRet;
    dLambda=uRet+3;
    du=uRet+6;
    c1=c;
    c2=c+1;
    c3=c+2;
    
    dBeta[0]=-rVal*sinBeta*cosLambda;
    dBeta[1]=-rVal*sinBeta*sinLambda;
    dBeta[2]=u*cosBeta;
    
    *c1=sqrt(dBeta[0]*dBeta[0]+dBeta[1]*dBeta[1]+dBeta[2]*dBeta[2]);
    for(i=0;i<3;i++){dBeta[i]=dBeta[i]/(*c1);}
      
    du[0]=u*cosBeta*cosLambda/rVal;
    du[1]=u*cosBeta*sinLambda/rVal;
    du[2]=sinBeta;
    
    *c3=sqrt(du[0]*du[0]+du[1]*du[1]+du[2]*du[2]);
    for(i=0;i<3;i++){du[i]=du[i]/(*c3);}

    dLambda[0]=-rVal*cosBeta*sinLambda;
    dLambda[1]=rVal*cosBeta*cosLambda;
    dLambda[2]=0.0;
    
    *c2=sqrt(dLambda[0]*dLambda[0]+dLambda[1]*dLambda[1]+dLambda[2]*dLambda[2]);

//If the point is too close to the poles, then it is possible that c2 is
//nearly equal to zero. However, a normalized dLambda can just be found by
//orthogonality: it is orthogonal to dBeta and du. dLambda=cross(dBeta,du);
    dLambda[0]=dBeta[1]*du[2]-dBeta[2]*du[1];
    dLambda[1]=dBeta[2]*du[0]-dBeta[0]*du[2];
    dLambda[2]=dBeta[0]*du[1]-dBeta[1]*du[0];
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
