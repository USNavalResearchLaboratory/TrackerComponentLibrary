/**CART2ELLIPSE A C-code implementation of a function to convert Cartesian
*               coordinates to ellipsoidal coordinates.
*
*INPUTS: cartPoints A matrix of the points in ECEF Cartesian coordinates
*                   that are to be transformed into ellipsoidal
*                   coordinates. Each column of cartPoints is of the
*                   format [x;y;z].
*           a       The semi-major axis of the reference ellipsoid. If
*                   this argument is omitted, the value in
*                   Constants.WGS84SemiMajorAxis is used.
*           f       The flattening factor of the reference ellipsoid. If
*                   this argument is omitted, the value in
*                   Constants.WGS84Flattening is used.
*
*OUTPUTS:   points  A matrix of the converted points. Each column of the
*                   matrix has the format [latitude;longitude;altitude],
*                   with latitude and longitude given in radians.
*
*The primary algorithm implemented in this function is [1], which is a
*modified version of the algorithm in [2].
*
*However, both techniques cited above will fail if the point in question
*is too close to the origin (deep within the Earth). In such an instance,
*Halley's method as described in [3] is used instead, with a minor
*modification. If one lets the algorithm run for an arbitrary number of
*iterations, there will generally be underflows, since the ratio of S and C
*matter, but both terms can drift by a constant factor during the
*iterations. Thus, after each iterative step, the values are normalized so
*that C=1 (it takes one division). Convergence of Fukushima's method is
*assumed to occur after six iterations.
*
* The algorithm can be compiled for use in Matlab  using the 
* CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*points=Cart2Ellipse(cartPoints,a,f);
*or
*points=Cart2Ellipse(cartPoints);
*
*REFERENCES:
*[1] I. Sofair "Improved method for calculating exact geodetic latitude and
*    altitude revisited," Journal of Guidance, Control, and Dynamics, vol.
*    23, no. 2, p. 369, Mar. 2000.
*[2] I. Sofair, "Improved method for calculating exact geodetic latitude
*    and altitude," Journal of Guidance, Control, and Dynamics, vol. 20,
*    no. 4, pp. 824-826, Jul.-Aug. 1997.
*[3] Fukushima, T., "Transformation from Cartesian to geodetic coordinates
*    accelerated by Halley's method", J.Geodesy (2006) 79: 689-693.
*
*April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*/
/*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.*/

/*This header is required by Matlab.*/
#include "mex.h"
//Needed for sqrt, fabs, isfinite, fmax, fmin, copysign, and atan2
#include <math.h>
#include "MexValidation.h"

static const double pi=3.1415926535897932384626433832795028841971693993751;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *points,a,b,b2,f,e2,eps2,ec;
    size_t numVec,i;
    mxArray *retMat;
    double *retData;
    
    if(nrhs>3||nrhs<1){
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
    
    points=(double*)mxGetData(prhs[0]);
    //points[0] is x
    //points[1] is y
    //points[2] is z
    
    if(nrhs>1&&!mxIsEmpty(prhs[1])) {
        a=getDoubleFromMatlab(prhs[1]);
    } else {
        a=getScalarMatlabClassConst("Constants", "WGS84SemiMajorAxis");
    }

    if(nrhs>2&&!mxIsEmpty(prhs[2])) {
        f=getDoubleFromMatlab(prhs[2]);
    } else {
        f=getScalarMatlabClassConst("Constants", "WGS84Flattening");   
    }
    
    b=a*(1-f);//The semi-minor axis of the reference ellipsoid.
    b2=b*b;

    //The square of the first numerical eccentricity. 
    e2=2*f-f*f;
    //The square of the second numerical eccentricity.
    eps2=a*a/(b2)-1;

    //This value is used if Fukushima's method is chosen.
    ec=sqrt(1-e2);

    //Allocate space for the return variables.
    retMat=mxCreateDoubleMatrix(3,numVec,mxREAL);
    retData=(double*)mxGetData(retMat);
    
    for(i=0;i<numVec;i++) {
        double *phi, *lambda, *h;
        double x0,y0,z0;
        double r0,p,s,q;
        
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
        if(q>=0) {//Use Sofair's algorithm
            double u,v,P,Q,t,c,w,z,Ne,val;
            
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
        } else {//Use Fukushima's algorithm.
            double Cc,P,Z,S,C;
            //A gets a value within the loop. This initialization is just
            //to silence a warning when compiling with
            //-Wconditional-uninitialized
            double A=0;
            size_t curIter;
            const size_t maxIter=6;
            
            P=r0/a;
            Z=(ec/a)*fabs(z0);

            S=Z;
            C=ec*P;

            //Loop until convergence. Assume convergence in 6 iterations.
            for(curIter=0;curIter<maxIter;curIter++) {
                double B,F,D,SNew,CNew;
                
                A=sqrt(S*S+C*C);
                B=1.5*e2*S*C*C*((P*S-Z*C)*A-e2*S*C);
                F=P*A*A*A-e2*C*C*C;
                D=Z*A*A*A+e2*S*S*S;

                SNew=D*F-B*S;
                CNew=F*F-B*C;

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
            }
            Cc=ec*C;

            //If the point is along the z-axis, then SNew and CNew will
            //both be zero, leading to a non-finite result.
            if(!isFinite(S)) {
                *phi=copySign(pi/2,z0);
                *h=fabs(z0)-b;
            } else {
                *phi=copySign(atan(S/Cc),z0);
                *h=(r0*Cc+fabs(z0)*S-b*A)/sqrt(Cc*Cc+S*S);
            }
        }
    }

    plhs[0]=retMat;
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
