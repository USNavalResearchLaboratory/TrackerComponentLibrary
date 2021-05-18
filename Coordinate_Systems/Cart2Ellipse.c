/**CART2ELLIPSE Convert Cartesian coordinates to ellipsoidal (latitude,
*               longitude, and altitude) coordinates.
*
*INPUTS: cartPoints A matrix of the points in ECEF Cartesian coordinates
*                   that are to be transformed into ellipsoidal
*                   coordinates. Each column of cartPoints is of the
*                   format [x;y;z].
*         algorithm This specified the algorithm to use for the conversion.
*                   Note that none work at the origin. Possible values are:
*                   0 (The default if this parameter is omitted or an empty
*                     matrix is passed and f<0.01) Use the algorithm of
*                     Olson in [1]. 
*                   1 Use the Algorithm of Sofair in [2], which is a
*                     modification of [3]. This will not work close to the
*                     center of the Earth.
*                   2 (The default if this parameter is omitted or an empty
*                     matrix is passed and f>=0.01) Use the algorithm of
*                     Fukushima in [4]. This should work close to the
*                     center of the Earth.
*                 a The semi-major axis of the reference ellipsoid. If this
*                   argument is omitted, the value in
*                   Constants.WGS84SemiMajorAxis is used.
*                 f The flattening factor of the reference ellipsoid. If
*                   this argument is omitted, the value in
*                   Constants.WGS84Flattening is used.
*
*OUTPUTS:   points  A matrix of the converted points. Each column of the
*                   matrix has the format [latitude;longitude;altitude],
*                   with latitude and longitude given in radians.
*
*The algorithm of Olson in [1] appears to be the most precise non-iterative
*method available. The method of Sofair in [2] and [3] is also a
*non-iterative algorithm, but tends to have singificantly worse accuracy.
*Fukushima's algorithm in [4] is iterative and is implemented assuming
*convergence in six or fewer iterations. Its accuracy appears to be
*marginally better than [1], but it is slower.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*points=Cart2Ellipse(cartPoints,algorithm,a,f);
*or
*points=Cart2Ellipse(cartPoints);
*
*REFERENCES:
*[1] D. K. Olson, "Converting Earth-centered, Earth-fixed coordinates to
*   geodetic coordinates," IEEE Transactions on Aerospace and Electronic
*   Systems, vol. 32, no. 1, pp. 473-476, Jan. 1996.
*[2] I. Sofair "Improved method for calculating exact geodetic latitude and
*    altitude revisited," Journal of Guidance, Control, and Dynamics, vol.
*    23, no. 2, p. 369, Mar. 2000.
*[3] I. Sofair, "Improved method for calculating exact geodetic latitude
*    and altitude," Journal of Guidance, Control, and Dynamics, vol. 20,
*    no. 4, pp. 824-826, Jul.-Aug. 1997.
*[4] Fukushima, T., "Transformation from Cartesian to geodetic coordinates
*    accelerated by Halley's method", Journal of Geodesy, vol. 79, no. 12,
*    pp. 689-693, Mar. 2006.
*
*October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
//Needed for sqrt, fabs, isfinite, fmax, fmin, copysign, and atan2
#include <math.h>
#include "MexValidation.h"

static const double pi=3.1415926535897932384626433832795028841971693993751;

//Function prototypes
void OlsonAlg(double *points,double *retData,size_t numVec,double a,double f);
int SofairAlg(double *points,double *retData,size_t numVec,double a,double f);
void FukishimaAlg(double *points,double *retData,size_t numVec,double a,double f);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *points,a,f;
    size_t numVec;
    mxArray *retMat;
    double *retData;
    int algorithm;
    
    if(nrhs>4||nrhs<1){
        mexErrMsgTxt("Wrong number of inputs");
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    checkRealDoubleArray(prhs[0]);
    numVec = mxGetN(prhs[0]);
    
    if(mxGetM(prhs[0])!=3) {
        mexErrMsgTxt("The input vector has a bad dimensionality.");
    }
    
    points=mxGetDoubles(prhs[0]);
    //points[0] is x
    //points[1] is y
    //points[2] is z
    
    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        a=getDoubleFromMatlab(prhs[2]);
    } else {
        a=getScalarMatlabClassConst("Constants", "WGS84SemiMajorAxis");
    }

    if(nrhs>3&&!mxIsEmpty(prhs[3])) {
        f=getDoubleFromMatlab(prhs[3]);
    } else {
        f=getScalarMatlabClassConst("Constants", "WGS84Flattening");   
    }
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        algorithm=getIntFromMatlab(prhs[1]);
    } else {
        if(f<0.01) {
            algorithm=0;
        } else {
            algorithm=2;
        }
    }

    //Allocate space for the return variables.
    retMat=mxCreateDoubleMatrix(3,numVec,mxREAL);
    retData=mxGetDoubles(retMat);
    
    switch(algorithm) {
        case 0://Olson's algorithm
            OlsonAlg(points,retData,numVec,a,f);
            break;
        case 1://Sofair's algorithm
        {
            int retVal;
            
            retVal=SofairAlg(points,retData,numVec,a,f);
            if(retVal!=0) {
                mexErrMsgTxt("The point given is too close to the center of the Earth for the algorithm of Sofair.");
            }
        }
        break;
        case 2://Fukushima's algorithm
            FukishimaAlg(points,retData,numVec,a,f);
            break;
        default:
            mexErrMsgTxt("Unknown algorithm specified.");
    }

    plhs[0]=retMat;
}

int SofairAlg(double *points,double *retData,size_t numVec,double a,double f) {
    size_t i;
    double b, b2, e2, eps2;
    
    b=a*(1-f);//The semi-minor axis of the reference ellipsoid.
    b2=b*b;

    //The square of the first numerical eccentricity. 
    e2=2*f-f*f;
    //The square of the second numerical eccentricity.
    eps2=a*a/(b2)-1;

    for(i=0;i<numVec;i++) {
        double *phi, *lambda, *h;
        double x0,y0,z0;
        double r0,p,s,q;
        double u,v,P,Q,t,c,w,z,Ne,val;
        
        //Get the Cartesian point to convert.
        x0=points[3*i];
        y0=points[3*i+1];
        z0=points[3*i+2];
        
        //Get the addresses of where the converted components will go
        phi=retData+3*i;
        lambda=retData+3*i+1;
        h=retData+3*i+2;
        
        r0=sqrt(x0*x0+y0*y0);
        p=fabs(z0)/eps2;
        s=r0*r0/(e2*eps2);
        q=p*p-b2+s;
    
        *lambda=atan2(y0,x0);
        
        //If the point is too deep within the Earth for this algorithm to
        //work.
        if(q<0) {
            return 1;
        }
  
        u=p/sqrt(q);
        v=b2*u*u/q;
        P=27.0*v*s/q;
        Q=pow(sqrt(P+1.0)+sqrt(P),2.0/3.0);
        t=(1.0+Q+1/Q)/6.0;
        c=u*u-1+2*t;
        //This condition prevents finite precision problems due to
        //subtraction within the square root.
        c=fMax(c,0);
        c=sqrt(c);
        w=(c-u)/2.0;

    //The z coordinate of the closest point projected on the ellipsoid.
    //The fmax command deals with precision problems when the argument
    //is nearly zero. The problems arise due to the subtraction within
    //the square root.
        z=sqrt(t*t+v)-u*w-t/2.0-1.0/4.0;
        z=fMax(z,0);
        z=copySign(sqrt(q)*(w+sqrt(z)),z0);

        Ne=a*sqrt(1+eps2*z*z/b2);

        //The min and max terms deals with finite precision problems.
        val=fMin(z*(eps2+1)/Ne,1);
        val=fMax(val,-1.0);
        *phi=asin(val);
        *h=r0*cos(*phi)+z0*sin(*phi)-a*a/Ne;
    }
    
    return 0;
}

void FukishimaAlg(double *points,double *retData,size_t numVec,double a,double f) {
    size_t i;
    double b, e2, ec;
    
    b=a*(1-f);//The semi-minor axis of the reference ellipsoid.

    //The square of the first numerical eccentricity. 
    e2=2*f-f*f;
    
    ec=sqrt(1-e2);
    
    for(i=0;i<numVec;i++) {
        double *phi, *lambda, *h;
        double x0,y0,z0;
        double r0;
        double Cc,P,Z,S,C;
        double A=0;
        size_t curIter;
        const size_t maxIter=500;
        
        //Get the Cartesian point to convert.
        x0=points[3*i];
        y0=points[3*i+1];
        z0=points[3*i+2];
        
        //Get the addresses of where the converted components will go
        phi=retData+3*i;
        lambda=retData+3*i+1;
        h=retData+3*i+2;
    
        r0=sqrt(x0*x0+y0*y0);
    
        *lambda=atan2(y0,x0);

        P=r0/a;
        Z=(ec/a)*fabs(z0);

        S=Z;
        C=ec*P;

        //Loop until convergence. Assume convergence in 6 iterations.
        for(curIter=0;curIter<maxIter;curIter++) {
            double B,F,D,SNew,CNew,SOld,COld;
            
            A=sqrt(S*S+C*C);
            B=1.5*e2*S*C*C*((P*S-Z*C)*A-e2*S*C);
            F=P*A*A*A-e2*C*C*C;
            D=Z*A*A*A+e2*S*S*S;

            SNew=D*F-B*S;
            CNew=F*F-B*C;
            
            SOld=S;
            COld=C;

            SNew=SNew/CNew;
            if(!isFinite(SNew)) {
                S=SNew;
                C=1;
                A=sqrt(S*S+C*C);
                break;
            } else {
                S=SNew;
                C=1;
            }
            
            if(SNew==SOld&&C==COld) {
                break;
            }
        }
        Cc=ec*C;
        //If the point is along the z-axis, then SNew and CNew will
        //both be zero, leading to a non-finite result.
        if(!isFinite(S)) {
            double temp=1.0;
            temp=copySign(temp,z0);
            *phi=temp*pi/2;
            *h=fabs(z0)-b;
        } else {
            double temp=1.0;
            temp=copySign(temp,z0);
            *phi=temp*atan(S/Cc);
            *h=(r0*Cc+fabs(z0)*S-b*A)/sqrt(Cc*Cc+S*S);
        }
    }
}

void OlsonAlg(double *points,double *retData,size_t numVec,double a,double f) {
    size_t i;
    double e2, a1, a2, a3, a4, a5, a6;
    //The square of the eccentricity.
    e2=2*f-f*f;

    a1=a*e2;
    a2=a1*a1;
    a3=a1*e2/2;
    a4=(5/2)*a2;
    a5=a1+a3;
    a6=1-e2;
    
    for(i=0;i<numVec;i++) {
        double *phi, *lambda, *h;
        double x,y,z;
        double zp,w2,w,z2,r2,r,s2,c2,s,c,ss;
        double g,rg,rf,u,v,m,f,p;
        
        //Get the Cartesian point to convert.
        x=points[3*i];
        y=points[3*i+1];
        z=points[3*i+2];
        
        //Get the addresses of where the converted components will go
        phi=retData+3*i;
        lambda=retData+3*i+1;
        h=retData+3*i+2;
        
        zp=fabs(z);
        w2=x*x+y*y;
        w=sqrt(w2);
        z2=z*z;
        r2=w2+z2;
        //The algorithm will work with points deep close to the origin.
        //Thus, there is no need to have a test for r being too small as
        //is the case in [1].
        r=sqrt(r2);

        *lambda=atan2(y,x);
        s2=z2/r2;
        c2=w2/r2;
        u=a2/r;
        v=a3-a4/r;
        if(c2>0.3) {
            s=(zp/r)*(1+c2*(a1+u+s2*v)/r);
            *phi=asin(s);
            ss=s*s;
            c=sqrt(1-ss);
        } else {
            c=(w/r)*(1-s2*(a5-u-c2*v)/r);
            *phi=acos(c);
            ss=1-c*c;
            s=sqrt(ss);
        }

        g=1-e2*ss;
        rg=a/sqrt(g);
        rf=a6*rg;
        u=w-rg*c;
        v=zp-rf*s;
        f=c*u+s*v;
        m=c*v-s*u;
        p=m/(rf/g+f);
        *phi=*phi+p;
        *h=f+m*p/2;
        if(z<0) {
        	*phi=-*phi;
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
