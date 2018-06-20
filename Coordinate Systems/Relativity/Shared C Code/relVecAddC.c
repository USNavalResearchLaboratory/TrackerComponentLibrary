/*RELVECADDC A C-only version of a function to add two 3X1 velocity vectors
 *           in a manner consistent with special relativity theory.
 *
 *INPUTS c   The speed of light as a double.
 *       v   A pointer to a 3D vector (double array) of the velocity  of
 *           the observer with respect to the inertial reference coordinate
 *           system.
 *       u   A pointer to a 3D vector (double array) of the velocity  of
 *           the of the object with respect to the observer's coordinate
 *           system.
 *   vObjRef A pointer to a 3D vector (double array) in which the result of
 *           the velocity addition will be placed.
 *
 *The formulae for special relativistic velocity addition is derived in
 *Chapter 1.4 of
 *G. Ludyk, Einstein in Matrix Form: Exact Derivation of the Theory of
 *Special and General Relativity without Tensors. Heidelberg: Springer,
 *2013.
 *The magnitudes of vObsFrame and vObjInFrame must both be less than the
 *speed of light.
 *
 *March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#include "relFuncs.h"
//For sqrt
#include <math.h>
//For size_t
#include <stddef.h>

void relVecAddC(double c, double *v,double *u, double *vObjRef) {
    double vMag=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    double gamma=1/sqrt(1-(vMag*vMag)/(c*c));
    double uv=v[0]*u[0]+v[1]*u[1]+v[2]*u[2];
    double Denom=1+uv/(c*c);
    size_t i;
    
    for(i=0;i<3;i++) {
        vObjRef[i]=(v[i]+u[i]+(1/gamma-1)*(u[i]-(uv/(vMag*vMag))*v[i]))/Denom;   
    }
    //If vMag=0 then NaNs will appear. We are not using the isnan function
    //to detect the NaNs, because many versions of Windows do not support
    //this function. However, if a value x is a NaN then x!=x should be
    //true. Given NaNs the correct solution is u, because the frame is not
    //moving.
    if(vObjRef[0]!=vObjRef[0]||vObjRef[1]!=vObjRef[1]||vObjRef[2]!=vObjRef[2]) {
        for(i=0;i<3;i++) {
            vObjRef[i]=u[i];   
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
