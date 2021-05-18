function [alpha1,dist,alpha2]=indirectGeodeticProb(latLonStart,latLonEnd,a,f)
%%INDIRECTGEODETICPROB Solve the indirect geodetic problem. That is, given
%                      two points on an ellipsoidal Earth, find the
%                      initial bearing and distance one must travel to
%                      take the shortest (geodesic) path between the
%                      points.
%
%INPUTS: latLonStart The 2X1 initial point given in geodetic latitude and
%                  longitude in radians of the format [latitude;longitude].
%                  The latitude must be between -pi/2 and pi/2 radians and
%                  the longitude between -pi and pi radians.
%        latLonEnd The 2X1  final point for the geodesic path given in 
%                  geodetic latitude  and longitude in radians. latLonEnd
%                  has the same format at latLonStart.
%                a The semi-major axis of the reference ellipsoid (in
%                  meters). If this argument is omitted or an empty matrix
%                  is passed, the value in Constants.WGS84SemiMajorAxis is
%                  used.
%                f The flattening factor of the reference ellipsoid. If
%                  this argument is omitted or an empty matrix is passed,
%                  the value in Constants.WGS84Flattening is used.
%
%OUTPUTS: alpha1 The scalar forward azimuth at the starting point in
%                radians East of true North on the reference ellipsoid.
%                This is the initial heading one would travel to go
%                between latLonStart and latLonEnd.
%           dist The geodetic distance between the starting and stopping
%                points in meters.
%         alpha2 The forward azimuth at the ending point in radians East
%                of true North on the reference ellipsoid.
%
%This function implements the algorithm of [1] using certain corrections
%and changes in [2] and [3]. Additional changes are commented in the code.
%
%In [1], series terms up to order 6 are given for a number of variables. In
%this implementation, values to the tenth order are given. Note that values
%to the 30th order have been posted online at
%https://geographiclib.sourceforge.io/html/geodseries30.html
%However, there is generally no point in using higher orders unless one
%uses higher than double-precision floating point arithmetic. The series
%are only valid for small values of f, which is appropriate when dealing
%with planets such as the Earth that are approximated with ellipsoids
%having a low eccentricity.
%
%In some parts of the implementation, trigonometric functions are directly
%used rather than identities avoiding transcendental functions. This is
%because it was observed that in Matlab, after running the code a few
%times, the trigonometric functions were actually faster. This could be an
%idiosyncrasy of Matlab or it could due to the transcendental instructions
%in x86 processors (FPATAN,FPSIN,FPCOS,FPSINCOS). However, the alpha1 value
%being estimated is computed avoiding trigonometric functions until the end
%because it was observed to improve finite precision errors.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%[2] C. F. F. Karney. (2013, 31 Aug.) Addenda and errata for papers on
%    geodesics. [Online].
%    Available: http://geographiclib.sourceforge.net/geod-addenda.html
%[3] C. F. F. Karney. (2011, 7 Feb.) Geodesics on an ellipsoid of
%    revolution.
%    [Online]. Available: http://arxiv.org/pdf/1102.1215.pdf
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

maxIters=17;
tinyVal=sqrt(realmin);
tolPol=2^14*eps(pi/2);%Used to determine whether a point is at a pole.
tolEq=eps();%Used to determine whether points are almost equatorial.
lambda12ConvergeEps=8*eps();

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

%Equation 3 of [1].
e2=f*(2-f);
%Equation 4 of [1].
ep2=e2/(1-e2);
%Obtaining the semi-minor axis from the semi-major axis and flattening
%factor.
b=a*(1-f);

%Equation 2 in [1].
n=f/(2-f);

phi1=latLonStart(1);
phi2=latLonEnd(1);
lambda1=latLonStart(2);
lambda2=latLonEnd(2);

%This is a fix for the case where both points are on the equator (or are
%almost on the equator). If both are too close to the equator, we offset
%them by a bit so as to avoid numeric problems with distances.
if(abs(phi1)<tolEq&&abs(phi2)<tolEq)
    %The inequalities for the sign avoid the issue that the sign functions
    %returns 0 if the input is zero.
    phi1Sign=2*(phi1>=0)-1;
    phi2Sign=2*(phi2>=0)-1;
    phi1=phi1Sign*tolEq;
    phi2=phi2Sign*tolEq;
end

%It is required that phi1<=0, phi1<=phi2<=-phi1. If this is not the case,
%we have to transform the points.

%For phi2<=-phi1, phi1 has to be the one with the larger magnitude and phi1
%has to be negative.
if(abs(phi1)<abs(phi2))
    %Swap points to make phi1 the one with the larger magnitude.
    temp=phi1;
    phi1=phi2;
    phi2=temp;
    
    temp=lambda1;
    lambda1=lambda2;
    lambda2=temp;
    swappedPoints=true;
else
    swappedPoints=false; 
end

if(phi1>0)
    phi1=-phi1;
    phi2=-phi2;
    swappedLatSigns=true;
else
    swappedLatSigns=false;
end
    
%It is required that 0<=lambda12<=pi. If that is not the case, flip the
%signs of the lambda values.
lambda12=wrapRange(lambda2-lambda1,-pi,pi);
if(abs(lambda12)<lambda12ConvergeEps)
    %This heuristic helps with finite precision limitations that can cause
    %massive errors. It is not described in Karney's papers.
    lambda12=0;
end

if(lambda12<0)
    %Flip the sign of the lambda difference so that 0<=lambda12<=pi.
    %Lambda1 and lambda2 are not reused, so they do not need to be
    %adjusted.
    lambda12=-lambda12;
    
    flippedLambdaSign=true;
else
    flippedLambdaSign=false;
end

%Now, 0<=lambda12<=pi and phi1<=0 and phi1<=phi2<=-phi1.

%Equation 6 in [1] for beta1 and beta2. Trajectories starting or ending
%near the poles are very slightly offset (by setting a minimum to cosBeta1)
%to avoid numeric issues.
beta1=atan((1-f)*tan(phi1));
beta2=atan((1-f)*tan(phi2));
sinBeta1=sin(beta1);
cosBeta1=max(tinyVal,cos(beta1));
sinBeta2=sin(beta2);
cosBeta2=max(tinyVal,cos(beta2));

if(lambda12==pi||lambda12==0||phi1+pi/2<=tolPol)
    %Check for the meridional special case.
    %alpha1=lambda12
    sinAlpha1=sin(lambda12);
    cosAlpha1=cos(lambda12);
    
    %It has already converged; we just need to tell it to compute all the
    %necessary values before exiting.
    deltaLambda12=0;
else
    %Not meridional.

    %The initial heading estimate so that we can start Newton's method.
    [sinAlpha1,cosAlpha1]=getInitialEstimate(lambda12,cosBeta1,cosBeta2,sinBeta1,sinBeta2,f,e2,ep2,n);
    %It has not already converged.
    deltaLambda12=Inf;
end

%The sine and cosine values of alpha1 that bracket the known region.
%Initially, the lower bound is alpha=0 and the upper bound is alpha=pi.
sinAlpha1LB=tinyVal;
cosAlpha1LB=1;
sinAlpha1UB=tinyVal;
cosAlpha1UB=-1;

for curIter=1:maxIters
    %Equation 10 in [1] for alpha0. We don't directly compute it.
    %Rather, we just compute the sine and cosine values. alpha0 is from
    %0 to pi.
    %alpha0=atan2(sinAlpha1*cosBeta1,sqrt(cosAlpha1^2+(sinAlpha1*sinBeta1)^2));
    %This method of computing the sine and cosine values is suggested
    %in Section 6 of [3].
    num=sinAlpha1*cosBeta1;
    denom=sqrt(cosAlpha1^2+(sinAlpha1*sinBeta1)^2);
    mag=hypot(num,denom);
    sinAlpha0=num/mag;
    cosAlpha0=denom/mag;

    %Equation 11 in [1] for sigma1:
    sigma1=atan2(sinBeta1,cosAlpha1*cosBeta1);
    sinSigma1=sin(sigma1);
    cosSigma1=cos(sigma1);
    sin2Sigma1=sin(2*sigma1);
    cos2Sigma1=cos(2*sigma1);

    %Equation 45 in [1] for cosine of alpha2. The max function deals with
    %possible finite precision errors.
    cosAlpha2=sqrt(max(0,cosAlpha1^2*cosBeta1^2+(cosBeta2^2-cosBeta1^2)))/cosBeta2;

    %Equation 11 in [1] for sigma2:
    sigma2=atan2(sinBeta2,cosAlpha2*cosBeta2);
    sinSigma2=sin(sigma2);
    cosSigma2=cos(sigma2);
    sin2Sigma2=sin(2*sigma2);
    cos2Sigma2=cos(2*sigma2);

    %Equation 9 in [1].
    k2=ep2*cosAlpha0^2;

    %Equation 16 in [1].
    temp=sqrt(1+k2);
    epsilon=(temp-1)/(temp+1);

    %Equations 17 and 18. These will be needed for Equation 15.
    [A1,C1]=computeA1C1(epsilon);
    %Equation 15 in [1].
    I1Sigma1=A1*(sigma1+evalSinCosSeries(C1,cos2Sigma1,sin2Sigma1));
    %Equation 15 in [1].
    I1Sigma2=A1*(sigma2+evalSinCosSeries(C1,cos2Sigma2,sin2Sigma2));

    if(abs(deltaLambda12)<=lambda12ConvergeEps||curIter==maxIters)
        alpha1=atan2(sinAlpha1,cosAlpha1);
        %Equation 5 in [1].
        sinAlpha2=sinAlpha0/cosBeta2;
        
        alpha2=atan2(sinAlpha2,cosAlpha2);
        %Equation 7 in [1].
        s1=b*I1Sigma1;
        %Equation 7 in [1].
        s2=b*I1Sigma2;
        dist=s2-s1;
        break;
    end

    %Equations 42 and 43 in [1], nominally, but actually Equation 42 was
    %replaced by a series with a (1+epsilon) term in front as suggested in
    %[2] for improved accuracy. These will be needed for Equation 41.
    [A2,C2]=computeA2C2(epsilon);
    %Equation 41 in [1].
    I2Sigma1=A2*(sigma1+evalSinCosSeries(C2,cos2Sigma1,sin2Sigma1));
    %Equation 41 in [1].
    I2Sigma2=A2*(sigma2+evalSinCosSeries(C2,cos2Sigma2,sin2Sigma2));

    %Equation 12 in [1].
    omega1=atan2(sinAlpha0*sinSigma1,cosSigma1);
    %Equation 12 in [1].
    omega2=atan2(sinAlpha0*sinSigma2,cosSigma2);

    %Equation 40 in [1].
    JSigma1=I1Sigma1-I2Sigma1;
    %Equation 40 in [1].
    JSigma2=I1Sigma2-I2Sigma2;

    %Equations 24 and 25 in [1].
    [A3Coeffs,C3Coeffs]=computeA3C3Coeffs(n);
    [A3,C3]=computeA3C3(A3Coeffs,C3Coeffs,epsilon);
    %Equation 23 in [1].
    I3Sigma1=A3*(sigma1+evalSinCosSeries(C3,cos2Sigma1,sin2Sigma1));
    %Equation 23 in [1].
    I3Sigma2=A3*(sigma2+evalSinCosSeries(C3,cos2Sigma2,sin2Sigma2));

    %Equation 8 in [1].
    lambda1Cur=omega1-f*sinAlpha0*I3Sigma1;
    %Equation 8 in [1].
    lambda2Cur=omega2-f*sinAlpha0*I3Sigma2;
    lambda12Cur=lambda2Cur-lambda1Cur;
    
    deltaLambda12=lambda12-lambda12Cur;

    %Equation 38 in [1].
    m12=b*(sqrt(1+k2*(sinSigma2)^2)*(cosSigma1*sinSigma2)-sqrt(1+k2*(sinSigma1)^2)*(sinSigma1*cosSigma2)-cosSigma1*cosSigma2*(JSigma2-JSigma1));
    dLambda12dAlpha1=(m12/a)*(1/(cosAlpha2*cosBeta2));

    if(~isfinite(dLambda12dAlpha1))
        %This should be the case where beta2=+/-beta1 and alpha1=pi/2.
        %The max is just in case of any finite precision issues.
        dLambda12dAlpha1=-2*(sqrt(max(0,1-e2*cosBeta1^2))/sinBeta1);
    end

    deltaAlpha1=deltaLambda12/dLambda12dAlpha1;
    %If we were directly updating the angle, the update would be:
    %alpha1=alpha1+deltaAlpha1;

    %Use angle sum formulae to update sunAlpha1 and cosAlpha1
    sinDeltaAlpha1=sin(deltaAlpha1);
    cosDeltaAlpha1=cos(deltaAlpha1);
    sinAlpha1New=sinAlpha1*cosDeltaAlpha1+cosAlpha1*sinDeltaAlpha1;
    cosAlpha1New=cosAlpha1*cosDeltaAlpha1-sinAlpha1*sinDeltaAlpha1;

    tanAlpha1New=cosAlpha1New/sinAlpha1New;
    if(deltaLambda12>0)
        %Update the upper bounds if the tangent (and hence the angle)
        %is larger than what has been seen thus far.
        tanAlpha1UB=cosAlpha1UB/sinAlpha1UB;
        if(tanAlpha1New>tanAlpha1UB)
            cosAlpha1UB=cosAlpha1New;
            sinAlpha1UB=sinAlpha1New;
        end
    elseif(deltaLambda12<0)
        %Update the lower bounds if the tangent (and hence the angle)
        %is less than what has been seen thus far.
        tanAlpha1LB=cosAlpha1LB/sinAlpha1LB;
        if(tanAlpha1New<tanAlpha1LB)
            cosAlpha1LB=cosAlpha1New;
            sinAlpha1LB=sinAlpha1New;
        end
    end

    %dLambda12dAlpha1>0 is equivalent to the "Hessian" in Newton's
    %method being positive definite (a descent direction). The
    %condition on sinAlpha1New is that it is bounded by 0 to pi. The
    %third condition is that the step is not too large.
    useAlternateUpdate=~(dLambda12dAlpha1>0&&sinAlpha1New>0&&abs(deltaAlpha1)<pi);

    %This test is described in [2] as a way of making the algorithm more
    %robust to poor initializations. Here, we have been keeping track
    %of valid upper bound and lower bound values of alpha (or better
    %put, of the sine and cosine of alpha). If the current update is
    %not valid, we will just average the current upper and lower
    %bounds (bisect the region).
    if(useAlternateUpdate)
        sinAlpha1=(sinAlpha1LB+sinAlpha1UB)/2;
        cosAlpha1=(cosAlpha1LB+cosAlpha1UB)/2;
    else
        sinAlpha1=sinAlpha1New;
        cosAlpha1=cosAlpha1New;
    end

    %Make sure the sine and cosine values are normalized.
    r=hypot(cosAlpha1,sinAlpha1);
    sinAlpha1=sinAlpha1/r;
    cosAlpha1=cosAlpha1/r;
end

%Finally, we have to undo the transformations.
if(swappedLatSigns)
    alpha1=-wrapRange(alpha1+pi,-pi,pi);
    alpha2=-wrapRange(alpha2+pi,-pi,pi);
end

if(flippedLambdaSign)
    alpha1=-alpha1;
    alpha2=-alpha2;
end

if(swappedPoints)
    temp=alpha1;
    alpha1=wrapRange(alpha2+pi,-pi,pi);
    alpha2=wrapRange(temp+pi,-pi,pi);
end
end

function [sinAlpha1,cosAlpha1]=getInitialEstimate(lambda12,cosBeta1,cosBeta2,sinBeta1,sinBeta2,f,e2,ep2,n)
%%GETINITIALESTIMATE Get the intiial estimate for Newton's method.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%[2] C. F. F. Karney. (2013, 31 Aug.) Addenda and errata for papers on
%    geodesics. [Online].
%    Available: http://geographiclib.sourceforge.net/geod-addenda.html
%[3] C. F. F. Karney. (2011, 7 Feb.) Geodesics on an ellipsoid of
%    revolution.
%    [Online]. Available: http://arxiv.org/pdf/1102.1215.pdf
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

%Tolerances for determining whether a point is nearly antipodal.
tolY=200*eps();
tolX=1000*sqrt(eps());

%As in the text before Equation 63 in [3], we shall rescale the coordinates
%first, to determine if we have a nearly antipodal scenario, and if it is
%a nearly antipodal scenario, to solve it better.

%As described in [3] before Equation 64, we need to get A3 assuming that
%alpha0=pi/2-abs(beta1). alpha0 comes into the calculation in the value of
%k in Equation 9 in [1]. Specifically, we need cos(alpha0). Given sinBeta1,
%we have cos(alpha0)=cos(pi/2-abs(asin(sinBeta1)))=abs(sinBeta1).
%Equation 9 in [1] (squared) with the above substitution.
k2=ep2*sinBeta1^2;

%Equation 16 in [1].
temp=sqrt(1+k2);
epsilon=(temp-1)/(temp+1);

%Equations 24 and 25 in [1].
A3Coeffs=computeA3C3Coeffs(n);
A3=computeA3C3(A3Coeffs,[],epsilon);

%From before Equation 64 in [1]:
deltaLambda=f*pi*A3*cosBeta1;
deltaBeta=cosBeta1*deltaLambda;
x=sin(lambda12-pi)/deltaLambda;

%Use an angle sum formula to get sin(beta1+beta2)
sinBeta1p2=sinBeta2*cosBeta1+cosBeta1* sinBeta1;
y=sinBeta1p2/deltaBeta;

%The antipodal point is now at (1,0). We are checking whether we are close
%to it.
nearlyAntipodal=(abs(y)<tolY&&abs(x-1)<tolX);

if(~nearlyAntipodal)
    %Equation 48 in [1]. The max is just in case of any finite precision
    %issues.
    wBar=sqrt(max(0,1-e2*((cosBeta1+cosBeta2)/2)^2));
    if(wBar==0)
        omega12=lambda12;
    else
        omega12=lambda12/wBar;
    end

    %Correction if the scaling is obviously invalid.
    if(abs(omega12)>pi)
        omega12=sign(omega12)*pi;
    end
    
    cosOmega12=cos(omega12);
    sinOmega12=sin(omega12);

    z1Real=cosBeta1*sinBeta2-sinBeta1*cosBeta2*cosOmega12;
    z1Imag=cosBeta2*sinOmega12;
    
    %Equation 49 in [1], but we return the sine and cosine separately
    %rather than using:
    %alpha1=atan2(z1Imag,z1Real);
    r=hypot(z1Imag,z1Real);
    sinAlpha1=z1Imag/r;
    cosAlpha1=z1Real/r;
else
    %Solve the astroid problem for nearly antipodal points.
    %We are using the scaled x and y on page 9 of [3].

    %Solve Equation 55 in [1].
    mu=solveAstroidEq(x,y);
    if(y~=0)
        %Equation 56 in [1], but we return the sine and cosine separately
        %rather than using:
        %alpha1=atan2(-x/(1+mu),y/mu);
        num=-x/(1+mu);
        denom=y/mu;
        r=hypot(num,denom);
        sinAlpha1=num/r;
        cosAlpha1=denom/r;
    else
        %We are finding the sine and cosine of:
        %alpha1=atan2(-x,sqrt(max(0,1-x^2)));
        num=-x;
        denom=sqrt(max(0,1-x^2));
        r=hypot(num,denom);
        sinAlpha1=num/r;
        cosAlpha1=denom/r;
    end
end
end

function k=solveAstroidEq(x,y)
%%SOLVEASTROIDEQ This function solves the equation
%                k^4+2*k^3-(p+q-1)*k^2-2*q*k-q=0 for k.
%                where p=x^2 and q=y^2 with x>=0 and y>=0.
%
%The technique of Appendix B of [1] using e=1 is implemented here.
%
%REFERENCES:
%[1] C. F. F. Karney. (2011, 7 Feb.) Geodesics on an ellipsoid of
%    revolution.
%    [Online]. Available: http://arxiv.org/pdf/1102.1215.pdf
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

p=x^2;
q=y^2;

r=(p+q-1)/6;

if(q~=0)
    S=p*q/4;
    d=S*(S+2*r^3);

    SR3Sum=S+r^3;
    sumSign=2*(SR3Sum>=0)-1;
    if(d>=0)
        T=nthroot(SR3Sum+sumSign*sqrt(d),3);

        if(T==0)
            u=0;
        else
            u=r+T+r^2/T;
        end
    else
        %T is complex, but u should be real.
        Psi=atan2(sqrt(-d), -SR3Sum);
        u=r*(1+2*cos(Psi/3));
    end
    v=sqrt(u^2+q);
    if(u<0)
        upmv=q/(v-u);%u+/-v
    else
        upmv=u+v;
    end
    %An abs is thrown in just in case of negativity, but it should not be
    %needed.
    w=abs((upmv-q)/(2*v));
    %k should be non-negative. upmv
    k=abs(upmv/(sqrt(upmv+w^2)+w));
else
    %If q=0, the problem reduces to k*(1+k)^2=k*p. If p=x^2, then the
    %solutions possible are 0, -1-x and -1+x. We shall just choose the
    %largest positive solution.
    
    vals=[0;-1-x;-1+x];
    k=max(vals);
end
end

function [A2,C2]=computeA2C2(e)
%%COMPUETA2C2 This function computes A2 and C2 as in Equation 17 of [1],
%            but the series has been extended to the tenth order. As noted
%            in [2], it is better to use the expansion for A2 with the
%            first term being (1+epsilon)^(-1) rather than
%            (1-epsilon)^(-1).
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%[2] C. F. F. Karney. (2013, 31 Aug.) Addenda and errata for papers on
%    geodesics. [Online].
%    Available: http://geographiclib.sourceforge.net/geod-addenda.html
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

A2Coeffs=[1;
          -3/4;
          -7/64;
          -11/256;
          -375/16384;
          -931/65536];

e2=e*e;
A2=polyValGen(A2Coeffs,e2,1)/(1+e);

C2Coeffs=[1/2,        3/16,       5/48,       35/512,     63/1280,    77/2048,        429/14336,  6435/262144,    12155/589824,   46189/2621440;
          1/16,       1/32,       5/256,      7/512,      21/2048,    33/4096,        429/65536,  715/131072,     0,              0;
          1/32,       35/2048,    23/2048,    133/16384,  51/8192,    2607/524288,    0,          0,              0,              0;
          41/2048,    47/4096,    191/24576,  47/8192,    0,          0,              0,          0,              0,              0;
          59/4096,    557/65536,  0,          0,          0,          0,              0,          0,              0,              0];
    
C2=zeros(10,1);
ePow=e;
C2(1)=ePow*polyValGen(C2Coeffs(:,1),e2,1);
ePow=e2;
C2(2)=ePow*polyValGen(C2Coeffs(:,2),e2,1);
ePow=ePow*e;
C2(3)=ePow*polyValGen(C2Coeffs(1:4,3),e2,1);
ePow=ePow*e;
C2(4)=ePow*polyValGen(C2Coeffs(1:4,4),e2,1);
ePow=ePow*e;
C2(5)=ePow*polyValGen(C2Coeffs(1:3,5),e2,1);
ePow=ePow*e;
C2(6)=ePow*polyValGen(C2Coeffs(1:3,6),e2,1);
ePow=ePow*e;
C2(7)=ePow*polyValGen(C2Coeffs(1:2,7),e2,1);
ePow=ePow*e;
C2(8)=ePow*polyValGen(C2Coeffs(1:2,8),e2,1);
ePow=ePow*e;
C2(9)=ePow*polyValGen(C2Coeffs(1,9),e2,1);
ePow=ePow*e;
C2(10)=ePow*polyValGen(C2Coeffs(1,10),e2,1);

end


function [A1,C1]=computeA1C1(e)
%COMPUTEA1C1 This function computes A1 and C1 as in Equation 17 of [1], but
%            the series has been extended to the tenth order.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Coefficients for A1 as in Equation 17, but extended to the tenth order.
A1Coeffs=[1;
          1/4;
          1/64;
          1/256;
          25/16384;
          49/65536];
e2=e*e;
A1=polyValGen(A1Coeffs,e2,1)/(1-e);

%Coefficients for the C1 terms as in Equation 17, but extended to the tenth
%order.
C1Coeffs=[-1/2,    -1/16,   -1/48,     -5/512,    -7/1280, -7/2048,     -33/14336, -429/262144, -715/589824, -2431/2621440;
           3/16,    1/32,    3/256,     3/512,     7/2048,  9/4096,      99/65536,  143/131072,  0,           0;
          -1/32,   -9/2048, -3/2048,   -11/16384, -3/8192, -117/524288,  0,         0,           0,           0;
           19/2048, 7/4096,  17/24576,  3/8192,    0,       0,           0,         0,           0,           0;
          -3/4096,  1/6553,  0,         0,         0,       0,           0,         0,           0,           0];
      
C1=zeros(10,1);
ePow=e;
C1(1)=ePow*polyValGen(C1Coeffs(:,1),e2,1);
ePow=e2;
C1(2)=ePow*polyValGen(C1Coeffs(:,2),e2,1);
ePow=ePow*e;
C1(3)=ePow*polyValGen(C1Coeffs(1:4,3),e2,1);
ePow=ePow*e;
C1(4)=ePow*polyValGen(C1Coeffs(1:4,4),e2,1);
ePow=ePow*e;
C1(5)=ePow*polyValGen(C1Coeffs(1:3,5),e2,1);
ePow=ePow*e;
C1(6)=ePow*polyValGen(C1Coeffs(1:3,6),e2,1);
ePow=ePow*e;
C1(7)=ePow*polyValGen(C1Coeffs(1:2,7),e2,1);
ePow=ePow*e;
C1(8)=ePow*polyValGen(C1Coeffs(1:2,8),e2,1);
ePow=ePow*e;
C1(9)=ePow*polyValGen(C1Coeffs(1,9),e2,1);
ePow=ePow*e;
C1(10)=ePow*polyValGen(C1Coeffs(1,10),e2,1);
end


function [A3Coeffs,C3Coeffs]=computeA3C3Coeffs(n)
%COMPUTEA3C3Coeffs This function computes the coefficients of the
%            polynomials in terms of epsilon for A3 and C3 as in Equations
%            24 and 25 of [1]. The series has been extended to the tenth
%            order.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

A3nCoeffs=-[-1,         0,           0,            0,          0,          0,            0,            0,           0,             0,            0;
            1/2,       -1/2,         0,            0,          0,          0,            0,            0,           0,             0,            0;
            1/4,        1/8,        -3/8,          0,          0,          0,            0,            0,           0,             0,            0;
            1/16,       3/16,        1/16,        -5/16,       0,          0,            0,            0,           0,             0,            0;
            3/64,       1/32,        5/32,         5/128,     -35/128,     0,            0,            0,           0,             0,            0;
            3/128,      5/128,       5/256,        35/256,     7/256,     -63/256,       0,            0,           0,             0,            0;
            5/256,      15/1024,     35/1024,      7/512,      63/512,     21/1024,     -231/1024,     0,           0,             0,            0;
            25/2048,    35/2048,     21/2048,      63/2048,    21/2048,    231/2048,     33/2048,     -429/2048,    0,             0,            0;
            175/16384,  35/4096,     63/4096,      63/8192,    231/8192,   33/4096,      429/4096,     429/32768,  -6435/32768,    0,            0;
            245/32768,  315/32768,   105/16384,    231/16384,  99/16384,   429/16384,    429/65536,    6435/65536,  715/65536,    -12155/65536,  0;
            441/65536,  735/131072,  1155/131072,  165/32768,  429/32768,  1287/262144,  6435/262144,  715/131072,  12155/131072,  2431/262144, -46189/262144]';

%Equation 24
A3Coeffs=zeros(11,1);
for k=1:11
    A3Coeffs(k)=polyValGen(A3nCoeffs(1:k,k),n,1);
end

if(nargout>1)
    C3Coeffs=zeros(11,10);
    C31nCoeffs=[1/4,         -1/4,         0,            0,          0,           0,          0,           0,           0,            0,            0;
               1/8,          0,          -1/8,          0,          0,           0,          0,           0,           0,            0,            0;
               3/64,         3/64,       -1/64,        -5/64,       0,           0,          0,           0,           0,            0,            0;
               5/128,        1/64,        1/64,        -1/64,      -7/128,       0,          0,           0,           0,            0,            0;
               3/128,        11/512,      3/512,        1/256,     -7/512,      -21/512,     0,           0,           0,            0,            0;
               21/1024,      5/512,       13/1024,      1/512,     -1/1024,     -3/256,     -33/1024,     0,           0,            0,            0;
               243/16384,    189/16384,   83/16384,     127/16384,  3/16384,    -51/16384,  -165/16384,  -429/16384,   0,            0,            0;
               435/32768,    109/16384,   1/128,        45/16384,   39/8192,    -11/16384,  -33/8192,    -143/16384,  -715/32768,    0,            0;
               345/32768,    953/131072,  259/65536,    365/65536,  95/65536,    47/16384,  -143/131072, -143/32768,  -1001/131072, -2431/131072,  0;
               2511/262144,  317/65536,   1355/262144,  165/65536,  531/131072,  89/131072,  107/65536,  -169/131072, -1157/262144, -221/32768,   -4199/262144]';
    for k=1:10
        C3Coeffs(k,1)=polyValGen(C31nCoeffs(1:(k+1),k),n,1);
    end

    C32nCoeffs=[1/16,        -3/32,         1/32,        0,           0,             0,             0,            0,            0,           0,            0;
               3/64,        -1/32,        -3/64,        1/32,        0,             0,             0,            0,            0,           0,            0;
               3/128,        1/128,       -9/256,      -3/128,       7/256,         0,             0,            0,            0,           0,            0;
               5/256,        1/256,       -1/128,      -7/256,      -3/256,         3/128,         0,            0,            0,           0,            0;
               27/2048,      69/8192,     -39/8192,    -47/4096,    -41/2048,      -45/8192,       165/8192,     0,            0,           0,            0;
               187/16384,    39/8192,      31/16384,   -63/8192,    -185/16384,    -119/8192,     -33/16384,     143/8192,     0,           0,            0;
               287/32768,    47/8192,      31/65536,   -3/2048,     -537/65536,    -41/4096,      -693/65536,    0,            1001/65536,  0,            0;
               255/32768,    249/65536,    43/16384,   -119/65536,  -25/8192,      -507/65536,    -35/4096,     -507/65536,    39/32768,    221/16384,    0;
               1675/262144,  2127/524288,  753/524288,  357/524288, -3109/1048576, -3873/1048576, -1821/262144, -1885/262144, -2977/524288, 1989/1048576, 12597/1048576]';
    for k=1:9
        C3Coeffs(k,2)=polyValGen(C32nCoeffs(1:(k+2),k),n,1);
    end

    C33nCoeffs=[5/192,       -3/64,        5/192,    -1/192,        0,            0,            0,            0,            0,            0,            0;
               3/128,       -5/192,      -1/64,      5/192,       -1/128,        0,            0,            0,            0,            0,            0;
               7/512,       -1/384,      -77/3072,   5/3072,       65/3072,     -9/1024,       0,            0,            0,            0,            0;
               3/256,       -1/1024,     -71/6144,  -47/3072,      9/1024,       25/1536,     -55/6144,      0,            0,            0,            0;
               139/16384,    143/49152,  -383/49152, -179/16384,  -121/16384,    547/49152,    605/49152,   -143/16384,    0,            0,            0;
               243/32768,    95/49152,   -41/16384,  -147/16384,  -389/49152,   -109/49152,    557/49152,    455/49152,   -273/32768,    0,            0;
               581/98304,    377/131072, -33/16384,  -907/196608, -515/65536,   -1937/393216,  89/98304,     2093/196608,  455/65536,   -1547/196608,  0;
               1383/262144,  103/49152,  -17/262144, -127/32768,  -3853/786432, -25/4096,     -2011/786432,  265/98304,    1895/196608,  85/16384,    -969/131072]';
    for k=1:8
        C3Coeffs(k,3)=polyValGen(C33nCoeffs(1:(k+3),k),n,1);
    end

    C34nCoeffs=[7/512,       -7/256,     5/256,       -7/1024,       1/1024,      0,            0,            0,            0,            0,           0;
               7/512,       -5/256,    -7/2048,       9/512,       -21/2048,     1/512,        0,            0,            0,            0,           0;
               9/1024,      -43/8192,  -129/8192,     39/4096,      91/8192,    -91/8192,      11/4096,      0,            0,            0,           0;
               127/16384,   -23/8192,  -165/16384,   -47/8192,      213/16384,   11/2048,     -175/16384,    13/4096,      0,            0,           0;
               193/32768,    3/8192,   -505/65536,   -227/32768,    75/65536,    801/65536,    165/131072,  -637/65536,    455/131072,   0,           0;
               171/32768,    25/65536, -259/65536,   -471/65536,   -351/131072,  605/131072,   41/4096,     -189/131072,  -1127/131072,  119/32768,   0;
               1121/262144, 339/262144, -801/262144, -2525/524288, -2519/524288, 73/131072,    1539/262144,  1989/262144, -1633/524288, -3927/524288, 969/262144]';
    for k=1:7
        C3Coeffs(k,4)=polyValGen(C34nCoeffs(1:(k+4),k),n,1);
    end

    C35nCoeffs=[21/2560,     -9/512,       15/1024,      -7/1024,    9/5120,       -1/5120,       0,             0,            0,           0,          0;
               9/1024,      -15/1024,     3/2048,        57/5120,  -5/512,         9/2560,      -1/2048,        0,            0,           0,          0;
               99/16384,    -91/16384,   -781/81920,     883/81920, 319/81920,    -783/81920,    387/81920,    -13/16384,     0,           0,          0;
               179/32768,   -55/16384,   -79/10240,     -27/81920,  461/40960,    -139/81920,   -65/8192,       441/81920,   -35/32768,    0,          0;
               141/32768,   -109/131072, -217/32768,    -219/65536, 1559/327680,   5431/655360, -203/40960,    -1943/327680,  369/65536,  -85/65536,   0;
               1013/262144, -15/32768,   -5399/1310720, -199/40960, 1267/1310720,  1007/163840,  6277/1310720, -527/81920,   -659/163840,  459/81920, -969/655360]';
    for k=1:6
        C3Coeffs(k,5)=polyValGen(C35nCoeffs(1:(k+5),k),n,1);
    end

    C36nCoeffs=[11/2048,     -99/8192,       275/24576,     -77/12288,     9/4096,        -11/24576,      1/24576,        0,             0,             0,            0;
               99/16384,    -275/24576,     55/16384,       167/24576,   -407/49152,      35/8192,      -55/49152,       1/8192,        0,             0,            0;
               143/32768,   -253/49152,    -1105/196608,    481/49152,   -73/196608,     -169/24576,     1067/196608,   -11/6144,       15/65536,      0,            0;
               33/8192,     -221/65536,    -23/4096,        457/196608,   267/32768,     -329/65536,    -69/16384,       375/65536,    -77/32768,      17/49152,     0;
               1711/524288, -4333/3145728, -16885/3145728, -1343/1572864, 17381/3145728,  8519/2097152, -42985/6291456, -4885/3145728,  8549/1572864, -5797/2097152, 969/2097152]';
    for k=1:5
        C3Coeffs(k,6)=polyValGen(C36nCoeffs(1:(k+6),k),n,1);
    end

    C37nCoeffs=[429/114688, -143/16384,  143/16384, -91/16384,   39/16384,    -11/16384,       13/114688,    -1/114688,       0,            0,           0;
               143/32768,  -143/16384,  65/16384,   65/16384,  -109/16384,    507/114688,    -27/16384,      39/114688,     -1/32768,      0,           0;
               429/131072, -299/65536, -13/4096,    269/32768, -601/229376,  -989/229376,     9475/1835008, -4667/1835008,  1157/1835008, -17/262144,   0;
               403/131072, -13/4096,   -521/131072, 393/114688, 1209/229376, -11001/1835008, -3979/3670016,  8821/1835008, -833/262144,    429/458752, -57/524288]';
    for k=1:4
        C3Coeffs(k,7)=polyValGen(C37nCoeffs(1:(k+7),k),n,1);
    end  

    C38nCoeffs=[715/262144, -429/65536,     455/65536,   -637/131072,   315/131072,   -55/65536,      13/65536,      -15/524288,    1/524288,      0,           0;
               429/131072, -455/65536,     1053/262144,  35/16384,    -1361/262144,   69/16384,     -2095/1048576,   77/131072,   -105/1048576,   1/131072,    0;
               663/262144, -4173/1048576, -1717/1048576, 3485/524288, -3825/1048576, -9469/4194304,  18469/4194304, -6137/2097152, 4455/4194304, -885/4194304, 19/1048576]';
    for k=1:3
        C3Coeffs(k,8)=polyValGen(C38nCoeffs(1:(k+8),k),n,1);
    end  

    C39nCoeffs=[2431/1179648, -663/131072,  1105/196608, -833/196608,  153/65536, -187/196608,   221/786432,   -15/262144,     17/2359296, -1/2359296,  0;
               663/262144,   -1105/196608, 1003/262144,  187/196608, -391/98304,  1003/262144, -3425/1572864,  1921/2359296, -13/65536,    17/589824, -1/524288]';
    for k=1:2
        C3Coeffs(k,9)=polyValGen(C39nCoeffs(1:(k+9),k),n,1);
    end

    C310nCoeffs=[4199/2621440, -4199/1048576, 4845/1048576, -969/262144, 2907/1310720, -10659/10485760, 741/2097152, -95/1048576, 17/1048576, -19/10485760, 1/10485760]';
    for k=1:1
        C3Coeffs(k,10)=polyValGen(C310nCoeffs(1:(k+10),k),n,1);
    end
end
end

function [A3,C3]=computeA3C3(A3Coeffs,C3Coeffs,e)
%COMPUTEA3C3 This function computes A3 and C3 as in Equations 24 and 25 of
%            [1], but the series has been extended to the tenth order. The
%            coefficients for A3 and C3 in terms of n must have been
%            already evaluated.
%
%REFERENCES:
%[1] C. F. F. Karney, "Algorithms for geodesics," Journal of Geodesy, vol.
%    87, no. 1, pp. 43-55, Jan. 2013.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

A3=polyValGen(A3Coeffs,e,1);

if(nargout>1)
    C3=zeros(10,1);
    ePow=e;
    C3(1)=ePow*polyValGen(C3Coeffs(1:10,1),e,1);
    ePow=ePow*e;
    C3(2)=ePow*polyValGen(C3Coeffs(1:9,2),e,1);
    ePow=ePow*e;
    C3(3)=ePow*polyValGen(C3Coeffs(1:8,3),e,1);
    ePow=ePow*e;
    C3(4)=ePow*polyValGen(C3Coeffs(1:7,4),e,1);
    ePow=ePow*e;
    C3(5)=ePow*polyValGen(C3Coeffs(1:6,5),e,1);
    ePow=ePow*e;
    C3(6)=ePow*polyValGen(C3Coeffs(1:5,6),e,1);
    ePow=ePow*e;
    C3(7)=ePow*polyValGen(C3Coeffs(1:4,7),e,1);
    ePow=ePow*e;
    C3(8)=ePow*polyValGen(C3Coeffs(1:3,8),e,1);
    ePow=ePow*e;
    C3(9)=ePow*polyValGen(C3Coeffs(1:2,9),e,1);
    ePow=ePow*e;
    C3(10)=ePow*polyValGen(C3Coeffs(1:1,10),e,1);
end
end

%LICENSE:
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
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
