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

void geogHeading2uVecENU2D(const size_t N, const double * const geoEastOfNorth, double *uENU);
void geogHeading2uVecENU3D(const size_t N, const double *const  geoEastOfNorth, const double * const angUpFromLevel, double *uENU);

void ruv2CartGenCPP(double *retData,const double *z,const bool useHalfRange,const double *zTx,const double *zRx,const double *M, bool hasW);
void Cart2RuvGenCPP(double *retData,const double *points,bool useHalfRange,double *zTx,double *zRx,double *M,bool includeW);

double getRangeRate2DCPP(const double *points,bool useHalfRange,const double *xTx,const double *xRx);
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

/*GETENUAXESCPP A C++ function to compute the basis vectors of an East-
 *              North-Up coordinate system.
 *
 *INPUTS: u, c Pointers to arrays in which the unit vectors for the ENU
 *             axes and the magnitudes of the unnormalized vectors are
 *             placed.
 *    plhPoint The 3X1 array at which the axes are to be found given in
 *             terms of [latitude;longitude;height] with the geodetic 
 *             latitude and longitude in radians and the height in
 *             meters. The latitude should be between -pi/2 and pi/2.
 * justVertical A boolean parameter. If true then u and c only for the Up
 *             direction will be returned.
 *           a The semi-major axis of the reference ellipsoid.
 *           f The flattening factor of the reference ellipsoid
 *
 *OUTPUTS: None. The results are placed in u and c.
 *
 *Further comments are given in the Matlab file getENUAxes.
 *
 *February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 **/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/
template<class T> 
void getENUAxesCPP(T *u, T *c,const T *plhPoint,const bool justVertical,const T a,const T f) {
    const T phi=plhPoint[0];//The latitude
    const T lambda=plhPoint[1];//The longitude
    const T h=plhPoint[2];//The height
    const T sinP=sin(phi);
    const T cosP=cos(phi);
    const T sinL=sin(lambda);
    const T cosL=cos(lambda);
    //Pointer to the normalized East, North, and Up vectors.
    T *u1,*u2,*u3;
    //Pointers to the magnitudes of the unnormalized ENU vectors.
    T *c1,*c2,*c3;
    size_t i;
    
    if(justVertical==false) {
        u1=u;
        u2=u+3;
        u3=u+6;
        c1=c;
        c2=c+1;
        c3=c+2;
    } else {
        u3=u;
        c3=c;
        //These variables are not used when justVertical==true, but are
        //initialized to get rid of a warning when compiled with 
        //-Wconditional-uninitialized
        u1=NULL;
        u2=NULL;
        c1=NULL;
        c2=NULL;
    }
    
    //u3 is dr/dh (Up)
    u3[0]=cosP*cosL;
    u3[1]=cosP*sinL;
    u3[2]=sinP;
    //Barring precision problems, this is always one.
    *c3=sqrt(u3[0]*u3[0]+u3[1]*u3[1]+u3[2]*u3[2]);
    for(i=0;i<3;i++){u3[i]/=(*c3);}

    if(justVertical==false) {
        //The square of the first numerical eccentricity
        const T e2=2*f-f*f;
        const T sqrtArg=sqrt(1-e2*sinP*sinP);
        //The normal radius of curvature.
        const T Ne=a/sqrtArg;
        //The derivative of the normal radius of curvature with respect to
        //phi.
        const T dNedPhi=a*e2*cosP*sinP/(sqrtArg*sqrtArg*sqrtArg);

        //u1 is dr/dlambda, normalized (East).
        u1[0]=-(Ne+h)*cosP*sinL;
        u1[1]=(Ne+h)*cosP*cosL;
        u1[2]=0;
        *c1=sqrt(u1[0]*u1[0]+u1[1]*u1[1]+u1[2]*u1[2]);

        //u2 is dr/dphi, normalized (North)
        const T temp=cosP*dNedPhi-(Ne+h)*sinP;
        u2[0]=temp*cosL;
        u2[1]=temp*sinL;
        u2[2]=(Ne*(1-e2)+h)*cosP+(1-e2)*dNedPhi*sinP;
        *c2=sqrt(u2[0]*u2[0]+u2[1]*u2[1]+u2[2]*u2[2]);
        for(i=0;i<3;i++){u2[i]=u2[i]/(*c2);}

        //If the point is too close to the poles, then it is possible that
        //c1 is nearly equal to zero. However, u1 can just be found by
        //orthogonality: it is orthogonal to u3 and u2. u1=cross(u2,u3);
        u1[0]=u2[1]*u3[2]-u2[2]*u3[1];
        u1[1]=u2[2]*u3[0]-u2[0]*u3[2];
        u1[2]=u2[0]*u3[1]-u2[1]*u3[0];
    }
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

template<class T>
void ECEF2ENUCPP(const size_t N, const T *tECEF, const T * const plhPoint,const T a, const T f,T *tENU, T *M) {
    T u[9];
    T c[3];
    T RxCart[3];
    
    getENUAxesCPP(u,c,plhPoint,false,a,f);
    
    //M is the transpose of u.
    for(size_t i1=0;i1<3;i1++) {
        for(size_t i2=0;i2<3;i2++) {
            M[i2+3*i1]=u[i1+3*i2];
        }
    }
    
	ellips2CartCPP(plhPoint[0],plhPoint[1],plhPoint[2], a, f, RxCart);

    for(size_t k=0;k<N;k++) {
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

template<class T>
void ENU2ECEFCPP(const size_t N, const T *tENU, const T * const plhPoint,const T a, const T f,T *tECEF, T *MInv) {
    T c[3];
    T RxCart[3];
    
    getENUAxesCPP(MInv,c,plhPoint,false,a,f);

	ellips2CartCPP(plhPoint[0],plhPoint[1],plhPoint[2], a, f, RxCart);

    for(size_t k=0;k<N;k++) {
        tECEF[0]=MInv[0]*tENU[0]+MInv[3]*tENU[1]+MInv[6]*tENU[2]+RxCart[0];
        tECEF[1]=MInv[1]*tENU[0]+MInv[4]*tENU[1]+MInv[7]*tENU[2]+RxCart[1];
        tECEF[2]=MInv[2]*tENU[0]+MInv[5]*tENU[1]+MInv[8]*tENU[2]+RxCart[2];
        
        tECEF+=3;
        tENU+=3;
    }
}

template <class T>
void geogHeading2uVecCPP(const T * const point, const size_t NGeo, const T *const geoEastOfNorth, const size_t NEl, const T * const angUpFromLevel, const T a, const T f, T *u) {
    T uLocal[9];
    T c[3];
    
    getENUAxesCPP(uLocal,c,point,false,a,f);
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
}

template <class T>
void uVec2GeogHeading(const T * const point,const size_t N, const T * u,const  T a,const T f, T *geoEastOfNorth, T *angUpFromLevel) {
    T uLocal[9];
    T c[3];
    
    getENUAxesCPP(uLocal,c,point,false,a,f);
    const T *uEast=uLocal;
    const T *uNorth=uLocal+3;
    const T *uUp=uLocal+6;
    
    //Perform dot products over all of the vectors.
    for(size_t k=0;k<N;k++) {
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
