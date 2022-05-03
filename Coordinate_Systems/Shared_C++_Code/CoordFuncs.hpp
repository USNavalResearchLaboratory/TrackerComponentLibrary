/**COORDFUNCS A header file for C++ implementations of coordinate
 *            conversion functions as well as Jacobians and Hessians
 *            related to coordinate conversions. See the files implementing
 *            each function for more details on their usage. Most of the
 *            documentation for the algorithms is with the corresponding
 *            Matlab implementations.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

#ifndef COORDFUNCS
#define COORDFUNCS
#include <cstddef>
//For sin and cos used in the templates below.
#include <cmath>

void spher2CartCPP(double *cartPoint,const double *point,size_t systemType);
void spher2CartGenCPP(double *retData,const double *point,size_t systemType,const bool useHalfRange,const double *zTx,const double *zRx,const double *M);
void Cart2SphereGenCPP(double *retData,const double *cartPoints,const size_t systemType,const bool useHalfRange,const double *zTx,const double *zRx,const double *M);
void spher2CartNoRangeCPP(double *retData,const double *point,const size_t systemType,const double *M);
void getEllipsHarmAxesCPP(double *u, double *c,const double *pointHarmon,const double E);
void Cart2EllipsHarmonCPP(double *pointsHarmon,const double *cartPoints, const size_t numPoints, const double E);

void ruv2CartGenCPP(double *retData,const double *z,const bool useHalfRange,const double *zTx,const double *zRx,const double *M, bool hasW);
void Cart2RuvGenCPP(double *retData,const double *points,bool useHalfRange,double *zTx,double *zRx,double *M,bool includeW);

double getRangeRate1DCPP(const double *xTar,bool useHalfRange,const double *xTx,const double *xRx);
double getRangeRate2DCPP(const double *xTar,bool useHalfRange,const double *xTx,const double *xRx);
double getRangeRate3DCPP(const double *xTar,bool useHalfRange,const double *xTx,const double *xRx);

void rangeGradientCPP(const size_t numRows,double *J,double *tempSpace,const double *point,const bool useHalfRange,const double *lTx,const double *lRx);
void spherAngGradientGenCPP(double *retMat,const double *xG,const size_t systemType,const double *lRx, const double *M);
void calcSpherConvJacobGenCPP(double *J,const double *zSpher,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M);
void calcSpherConvJacobCPP(double *J, const double *point, const size_t systemType);
void calcSpherInvJacobCPP(double *J,const double *z,const size_t systemType);

void rangeHessianCPP(const size_t numDim,double *H,const double *x,const bool useHalfRange);
void rangeHessianGenCPP(const size_t numDim,double *H,const double *x,const bool useHalfRange,const double *lTx,const double *lRx);
void spherAngHessianCPP(double *HTotal,const double *xG,const size_t systemType);
void spherAngHessianGenCPP(double *HTotal,const double *xG,const size_t systemType,const double *lRx,const double *M);
void calcSpherInvHessianCPP(double *HTotal,const double *z,const size_t systemType);
void spherAngHessianGenCPP(double *H,const double *x,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M);
void calcSpherHessianCPP(double *H,const double *x,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M);
void calcSpherConvHessianCPP(double *H, const double *point, const size_t systemType, const bool useHalfRange);
void calcSpherConvHessianGenCPP(double *H,const double *zSpher,const size_t systemType,const bool useHalfRange,const double *lTx,const double *lRx,const double *M);

/**GETENUAXISDIRECTIONS A C++ function to compute the unit basis vectors of
 *              an East-North-Up coordinate system.
 *
 *INPUTS: u A pointer to an array in which the unit vectors for the ENU
 *          axes are placed. The first three elements are the first vector,
 *          the next three the second one, etc. If justVertical=true, then
 *          only a single vector, the vertical is returned.
 * plhPoint A length 2 array at which the axes are to be found given in
 *          terms of [latitude;longitude] with the geodetic latitude and
 *          longitude in radians and the height in meters. The latitude
 *          should be between -pi/2 and pi/2.
 * justVertical A boolean parameter. If true then u and c only for the Up
 *          direction will be returned.
 *
 *March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T> 
void getENUAxisDirections(T *u ,const T *plhPoint,const bool justVertical) { 
    const T phi=plhPoint[0];//The latitude
    const T lambda=plhPoint[1];//The longitude
    const T sinP=sin(phi);
    const T cosP=cos(phi);
    const T sinL=sin(lambda);
    const T cosL=cos(lambda);
    
    if(justVertical==false) {
        //u1 is dr/dlambda, normalized (East).
        u[0]=-sinL;
        u[1]=cosL;
        u[2]=0.0;
        u+=3;

        //u2 is dr/dphi, normalized (North)
        u[0]=-cosL*sinP;
        u[1]=-sinL*sinP;
        u[2]=cosP;
        u+=3;
    }
    //u3 is dr/dh (Up)
    u[0]=cosP*cosL;
    u[1]=cosP*sinL;
    u[2]=sinP;
}

/**GETENUAXESDIRMAG A C++ function to compute the basis vectors of an East-
 *              North-Up coordinate system as well as the magnitudes of the
 *              unnormalized vectors. The unnormalized vectors are
 *              gradients with respect to ellipsoidal coordinates.
 *
 *INPUTS: u A pointer to an array in which the unit vectors for the ENU
 *          axes are placed. The first three elements are the first vector,
 *          the next three the second one, etc.
 *        c A pointer to a length-3 array in which the magnitudes of the
 *          unnormalized vectors are placed.
 * plhPoint A length 2 array at which the axes are to be found given in
 *          terms of [latitude;longitude] with the geodetic latitude and
 *          longitude in radians and the height in meters. The latitude
 *          should be between -pi/2 and pi/2.
 *        a The semi-major axis of the reference ellipsoid.
 *        f The flattening factor of the reference ellipsoid
 *
 *March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T> 
void getENUAxesDirMag(T *u, T *c, const T *plhPoint, const T a, const T f) {
    const T phi=plhPoint[0];//The latitude
    const T lambda=plhPoint[1];//The longitude
    const T h=plhPoint[2];//The ellipsoidal height.
    const T sinP=sin(phi);
    const T cosP=cos(phi);
    const T sinL=sin(lambda);
    const T cosL=cos(lambda);
    //The square of the first numerical eccentricity
    const T e2=2*f-f*f;
    const T sqrtArg=sqrt(1-e2*sinP*sinP);
    //The normal radius of curvature.
    const T Ne=a/sqrtArg;
    //The derivative of the normal radius of curvature with respect to
    //phi.
    const T dNedPhi=a*e2*cosP*sinP/(sqrtArg*sqrtArg*sqrtArg);
    
    T ca,cb,cc;
    const T temp=cosP*dNedPhi-(Ne+h)*sinP;
    //u1 is dr/dlambda, normalized (East).
    u[0]=-sinL;
    u[1]=cosL;
    u[2]=0.0;
    u+=3;

    //Get the mangitude of the East vector.
    ca=-(Ne+h)*cosP*sinL;
    cb=(Ne+h)*cosP*cosL;
    c[0]=sqrt(ca*ca+cb*cb);

    //u2 is dr/dphi, normalized (North)
    u[0]=-cosL*sinP;
    u[1]=-sinL*sinP;
    u[2]=cosP;
    u+=3;
        
    //Get the magnitude of the North vector.
    ca=temp*cosL;
    cb=temp*sinL;
    cc=(Ne*(1-e2)+h)*cosP+(1-e2)*dNedPhi*sinP;
    c[1]=sqrt(ca*ca+cb*cb+cc*cc);
    //u3 is dr/dh (Up)
    u[0]=cosP*cosL;
    u[1]=cosP*sinL;
    u[2]=sinP;
    
    c[2]=1.0;
}

template<class T>
void ellips2CartCPP(const T latitude, const T longitude, const T height, const T a, const T f, T *cartPoint) {

    //The square of the first numerical eccentricity
    const T e2=2*f-f*f;
    
    //The geodetic latitude
    const T phi=latitude;
    //The longitude.
    const T lambda=longitude;
    //The altitude.
    const T h=height;

    const T sinP=sin(phi);
    const T cosP=cos(phi);
    const T sinL=sin(lambda);
    const T cosL=cos(lambda);

    //The normal radius of curvature.
    const T Ne=a/sqrt(1-e2*sinP*sinP);

    cartPoint[0]=(Ne+h)*cosP*cosL;
    cartPoint[1]=(Ne+h)*cosP*sinL;
    cartPoint[2]=(Ne*(1-e2)+h)*sinP;
}

/**ECEF2ENUCPP Convert points from ECEF coordinates to ENU coordinates.
 *
 *March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
void ECEF2ENUCPP(const size_t numPts, const T *tECEF, const T * const plhPoint,const T a, const T f,T *tENU, T *M) {
    T u[9];
    T RxCart[3];
    
    getENUAxisDirections(u ,plhPoint,false); 
    
    //M is the transpose of u.
    for(size_t i1=0;i1<3;i1++) {
        for(size_t i2=0;i2<3;i2++) {
            M[i2+3*i1]=u[i1+3*i2];
        }
    }
    
	ellips2CartCPP(plhPoint[0],plhPoint[1],plhPoint[2], a, f, RxCart);

    for(size_t k=0;k<numPts;k++) {
        T deltaVals[3];

        deltaVals[0]=tECEF[0]-RxCart[0];
        deltaVals[1]=tECEF[1]-RxCart[1];
        deltaVals[2]=tECEF[2]-RxCart[2];
        
        tENU[0]=M[0]*deltaVals[0]+M[3]*deltaVals[1]+M[6]*deltaVals[2];
        tENU[1]=M[1]*deltaVals[0]+M[4]*deltaVals[1]+M[7]*deltaVals[2];
        tENU[2]=M[2]*deltaVals[0]+M[5]*deltaVals[1]+M[8]*deltaVals[2];
        
        tECEF+=3;
        tENU+=3;
    }
}

/**ENU2ECEF Convert points from ECU coordinates to ECEF coordinates.
 *
 *March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
template<class T>
void ENU2ECEFCPP(const size_t numPoints, const T *tENU, const T * const plhPoint,const T a, const T f,T *tECEF, T *MInv) {
    T RxCart[3];

    getENUAxisDirections(MInv, plhPoint,false); 

	ellips2CartCPP(plhPoint[0],plhPoint[1],plhPoint[2], a, f, RxCart);

    for(size_t k=0;k<numPoints;k++) {
        tECEF[0]=MInv[0]*tENU[0]+MInv[3]*tENU[1]+MInv[6]*tENU[2]+RxCart[0];
        tECEF[1]=MInv[1]*tENU[0]+MInv[4]*tENU[1]+MInv[7]*tENU[2]+RxCart[1];
        tECEF[2]=MInv[2]*tENU[0]+MInv[5]*tENU[1]+MInv[8]*tENU[2]+RxCart[2];
        
        tECEF+=3;
        tENU+=3;
    }
}

template <class T>
void geogHeading2uVecCPP(const size_t NPts, const T * const point, const size_t NGeo, const T *const geoEastOfNorth, const size_t NEl, const T * const angUpFromLevel, T *u) {
    T uLocal[9];
    
    if(NPts==1) {
        getENUAxisDirections(uLocal, point,false); 
        
        const T *uEast=uLocal;
        const T *uNorth=uLocal+3;
        const T *uUp=uLocal+6;
        
        if(NEl>1&&NGeo>1) {//Assume both equal NGEo
            for(size_t i=0;i<NGeo;i++) {
                const T cosEl=cos(angUpFromLevel[i]);
                const T sinEl=sin(angUpFromLevel[i]);
                const T sinGeo=sin(geoEastOfNorth[i]);
                const T cosGeo=cos(geoEastOfNorth[i]);
                const T prodE=sinGeo*cosEl;
                const T prodN=cosGeo*cosEl;
     
                u[0]=prodE*uEast[0]+prodN*uNorth[0]+sinEl*uUp[0];
                u[1]=prodE*uEast[1]+prodN*uNorth[1]+sinEl*uUp[1];
                u[2]=prodE*uEast[2]+prodN*uNorth[2]+sinEl*uUp[2];
                u+=3;
            }
        } else if(NEl>1) {//NEl>1 and NGeo=1
            const T sinGeo=sin(geoEastOfNorth[0]);
            const T cosGeo=cos(geoEastOfNorth[0]);
            
            for(size_t i=0;i<NEl;i++) {
                const T cosEl=cos(angUpFromLevel[i]);
                const T sinEl=sin(angUpFromLevel[i]);
                const T prodE=sinGeo*cosEl;
                const T prodN=cosGeo*cosEl;
     
                u[0]=prodE*uEast[0]+prodN*uNorth[0]+sinEl*uUp[0];
                u[1]=prodE*uEast[1]+prodN*uNorth[1]+sinEl*uUp[1];
                u[2]=prodE*uEast[2]+prodN*uNorth[2]+sinEl*uUp[2];
                u+=3;
            }
        } else {//NEl=1 and NGeo>=1
            const T cosEl=cos(angUpFromLevel[0]);
            const T sinEl=sin(angUpFromLevel[0]);
            for(size_t i=0;i<NGeo;i++) {
                const T sinGeo=sin(geoEastOfNorth[i]);
                const T cosGeo=cos(geoEastOfNorth[i]);
                const T prodE=sinGeo*cosEl;
                const T prodN=cosGeo*cosEl;
     
                u[0]=prodE*uEast[0]+prodN*uNorth[0]+sinEl*uUp[0];
                u[1]=prodE*uEast[1]+prodN*uNorth[1]+sinEl*uUp[1];
                u[2]=prodE*uEast[2]+prodN*uNorth[2]+sinEl*uUp[2];
                u+=3;
            }
        }
    } else {
        if(NEl>1&&NGeo>1) {//Assume both equal NPts
            for(size_t i=0;i<NPts;i++) {
                getENUAxisDirections(uLocal, point+3*i,false); 
                const T *uEast=uLocal;
                const T *uNorth=uLocal+3;
                const T *uUp=uLocal+6;
    
                const T cosEl=cos(angUpFromLevel[i]);
                const T sinEl=sin(angUpFromLevel[i]);
                const T sinGeo=sin(geoEastOfNorth[i]);
                const T cosGeo=cos(geoEastOfNorth[i]);
                const T prodE=sinGeo*cosEl;
                const T prodN=cosGeo*cosEl;
     
                u[0]=prodE*uEast[0]+prodN*uNorth[0]+sinEl*uUp[0];
                u[1]=prodE*uEast[1]+prodN*uNorth[1]+sinEl*uUp[1];
                u[2]=prodE*uEast[2]+prodN*uNorth[2]+sinEl*uUp[2];
                u+=3;
            }
        } else if(NEl>1) {//NEl>1 and NGeo=1
            const T sinGeo=sin(geoEastOfNorth[0]);
            const T cosGeo=cos(geoEastOfNorth[0]);

            for(size_t i=0;i<NPts;i++) {
                getENUAxisDirections(uLocal, point+3*i,false); 
                const T *uEast=uLocal;
                const T *uNorth=uLocal+3;
                const T *uUp=uLocal+6;
    
                const T cosEl=cos(angUpFromLevel[i]);
                const T sinEl=sin(angUpFromLevel[i]);
                const T prodE=sinGeo*cosEl;
                const T prodN=cosGeo*cosEl;
     
                u[0]=prodE*uEast[0]+prodN*uNorth[0]+sinEl*uUp[0];
                u[1]=prodE*uEast[1]+prodN*uNorth[1]+sinEl*uUp[1];
                u[2]=prodE*uEast[2]+prodN*uNorth[2]+sinEl*uUp[2];
                u+=3;
            }
        } else if(NGeo>1) {
            const T cosEl=cos(angUpFromLevel[0]);
            const T sinEl=sin(angUpFromLevel[0]);

            for(size_t i=0;i<NPts;i++) {
                getENUAxisDirections(uLocal, point+3*i,false); 
                const T *uEast=uLocal;
                const T *uNorth=uLocal+3;
                const T *uUp=uLocal+6;
 
                const T sinGeo=sin(geoEastOfNorth[i]);
                const T cosGeo=cos(geoEastOfNorth[i]);
                const T prodE=sinGeo*cosEl;
                const T prodN=cosGeo*cosEl;
     
                u[0]=prodE*uEast[0]+prodN*uNorth[0]+sinEl*uUp[0];
                u[1]=prodE*uEast[1]+prodN*uNorth[1]+sinEl*uUp[1];
                u[2]=prodE*uEast[2]+prodN*uNorth[2]+sinEl*uUp[2];
                u+=3;
            }
        } else {//NEl=1 and NGeo=1
            const T sinGeo=sin(geoEastOfNorth[0]);
            const T cosGeo=cos(geoEastOfNorth[0]);
            const T cosEl=cos(angUpFromLevel[0]);
            const T sinEl=sin(angUpFromLevel[0]);
            const T prodE=sinGeo*cosEl;
            const T prodN=cosGeo*cosEl;

            for(size_t i=0;i<NPts;i++) {
                getENUAxisDirections(uLocal, point+3*i,false); 
                const T *uEast=uLocal;
                const T *uNorth=uLocal+3;
                const T *uUp=uLocal+6;
 
                u[0]=prodE*uEast[0]+prodN*uNorth[0]+sinEl*uUp[0];
                u[1]=prodE*uEast[1]+prodN*uNorth[1]+sinEl*uUp[1];
                u[2]=prodE*uEast[2]+prodN*uNorth[2]+sinEl*uUp[2];
                u+=3;
            }
        }
    }
}

/**UVEC2GEOGHEADING Obtain the geographic headings in radians East of true
*                 North as well as the elevations above the local tangent
*                 plane according to a particular reference ellipsoid that
*                 correspond to directions specified by unit vectors in
*                 ECEF coordinates.
*
*March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
**/
template <class T>
void uVec2GeogHeading(const T * const point,const size_t numPoints, const T * u, T *geoEastOfNorth, T *angUpFromLevel) {
    T uLocal[9];
    
    getENUAxisDirections(uLocal, point,false); 
    const T *uEast=uLocal;
    const T *uNorth=uLocal+3;
    const T *uUp=uLocal+6;

    for(size_t k=0;k<numPoints;k++) {
        //These are dot products.
        const T magEast=u[0]*uEast[0]+u[1]*uEast[1]+u[2]*uEast[2];
        const T magNorth=u[0]*uNorth[0]+u[1]*uNorth[1]+u[2]*uNorth[2];
        const T magUp=u[0]*uUp[0]+u[1]*uUp[1]+u[2]*uUp[2];
 
        geoEastOfNorth[k]=atan2(magEast,magNorth);
        angUpFromLevel[k]=asin(magUp);
        
        u+=3;
    }
}

#endif

/*LICENSE:
*
*The source code is in the public domain and not licensed or under
*copyright. The information and software may be used freely by the public.
*As required by 17 U.S.C. 403, third parties producing copyrighted works
*consisting predominantly of the material produced by U.S. government
*agencies must provide notice with such work(s) identifying the U.S.
*Government material incorporated and stating that such material is not
*subject to copyright protection.
*
*Derived works shall not identify themselves in a manner that implies an
*endorsement by or an affiliation with the Naval Research Laboratory.
*
*RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
*SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
*RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
*OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
