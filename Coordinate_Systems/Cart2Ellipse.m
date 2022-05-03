function points=Cart2Ellipse(cartPoints,algorithm,a,f)
%%CART2ELLIPSE Convert Cartesian coordinates to ellipsoidal (latitude,
%              longitude, and altitude) coordinates.
%
%INPUTS: cartPoints A matrix of the points in ECEF Cartesian coordinates
%                   that are to be transformed into ellipsoidal
%                   coordinates. Each column of cartPoints is of the
%                   format [x;y;z].
%         algorithm This specified the algorithm to use for the conversion.
%                   Note that none work at the origin. Possible values are:
%                   0 (The default if this parameter is omitted or an empty
%                     matrix is passed and f<0.01) Use the algorithm of
%                     Olson in [1]. 
%                   1 Use the Algorithm of Sofair in [2], which is a
%                     modification of [3]. This will not work close to the
%                     center of the Earth.
%                   2 (The default if this parameter is omitted or an empty
%                     matrix is passed and f>=0.01) Use the algorithm of
%                     Fukushima in [4]. This should work close to the
%                     center of the Earth.
%                 a The semi-major axis of the reference ellipsoid. If this
%                   argument is omitted, the value in
%                   Constants.WGS84SemiMajorAxis is used.
%                 f The flattening factor of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84Flattening is used.
%
%OUTPUTS: points A matrix of the converted points. Each column of the
%                matrix has the format [latitude;longitude;altitude], with
%                latitude and longitude given in radians.
%
%The algorithm of Olson in [1] appears to be the most precise non-iterative
%method available for targets far above the Earth. The method of Sofair in
%[2] and [3] is also a non-iterative algorithm, but tends to have
%significantly worse accuracy for such targets. Fukushima's algorithm in
%[4] is iterative and typically converges in six or fewer iterations. It is
%set to run for a maximum number of 500 iterations. More than 6 iterations
%can be necessary when a large flattening is used and the point in question
%is near the center of the Earth. Its accuracy appears to be marginally
%better than [1] for things not near the surface of the Earth/ reference
%ellipsoid, but it is slower.
%
%REFERENCES:
%[1] D. K. Olson, "Converting Earth-centered, Earth-fixed coordinates to
%    geodetic coordinates," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 32, no. 1, pp. 473-476, Jan. 1996.
%[2] I. Sofair "Improved method for calculating exact geodetic latitude and
%    altitude revisited," Journal of Guidance, Control, and Dynamics, vol.
%    23, no. 2, p. 369, Mar. 2000.
%[3] I. Sofair, "Improved method for calculating exact geodetic latitude
%    and altitude," Journal of Guidance, Control, and Dynamics, vol. 20,
%    no. 4, pp. 824-826, Jul.-Aug. 1997.
%[4] Fukushima, T., "Transformation from Cartesian to geodetic coordinates
%    accelerated by Halley's method", Journal of Geodesy, vol. 79, no. 12,
%    pp. 689-693, Mar. 2006.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<2||isempty(algorithm))
    if(f<0.01)
        algorithm=0;
    else
        algorithm=2;
    end
end

switch(algorithm)
    case 0%Olson's algorithm
        [lambda,phi,h]=OlsonAlg(cartPoints,a,f);
        
        if(any(imag(lambda)~=0)||any(imag(phi)~=0)||any(imag(phi)~=0))
            error('The point given is too close to the center of the Earth for the algorithm of Olson.')
        end
        
    case 1%Sofair's algorithm
        [lambda,phi,h]=SofairAlg(cartPoints,a,f);
    case 2%Fukushima's algorithm
        [lambda,phi,h]=FukishimaAlg(cartPoints,a,f);
    otherwise
        error('Unknown algorithm specified.')
end

points=[phi;lambda;h];

end

function [lambda,phi,h]=SofairAlg(cartPoints,a,f)
%%SOFAIRALG This implements the algorithm of [1], which is a modified
%           version of the algorithm of [2]. Both techniques will fail if
%           the point in question is too close to the origin (deep within
%           the Earth). If the algorithm fails, then this function wil have
%           an error.
%
%REFERENCES:
%[1] I. Sofair "Improved method for calculating exact geodetic latitude and
%    altitude revisited," Journal of Guidance, Control, and Dynamics, vol.
%    23, no. 2, p. 369, Mar. 2000.
%[2] I. Sofair, "Improved method for calculating exact geodetic latitude
%    and altitude," Journal of Guidance, Control, and Dynamics, vol. 20,
%    no. 4, pp. 824-826, Jul.-Aug. 1997.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPoints=size(cartPoints,2);

%The semi-minor axis of the reference ellipsoid.
b=a*(1-f);

%The square of the first numerical eccentricity. 
e2=2*f-f^2;

%The square of the second numerical eccentricity.
eps2=a^2/b^2-1;

%Allocate space for the results.
phi=zeros(1,numPoints);
lambda=zeros(1,numPoints);
h=zeros(1,numPoints);
for curPoint=1:numPoints
    %Extract the coordinates
    x0=cartPoints(1,curPoint);
    y0=cartPoints(2,curPoint);
    z0=cartPoints(3,curPoint);
    
    r0=sqrt(x0.^2+y0.^2);
    p=abs(z0)/eps2;
    s=r0.^2/(e2*eps2);
    q=p.^2-b.^2+s;
    
    lambda(curPoint)=atan2(y0,x0);
    
    if(q<0)
        error('The point given is too close to the center of the Earth for the algorithm of Sofair.')
    end

    u=p./sqrt(q);
    v=b^2*u.^2./q;
    P=27*v.*s./q;
    Q=(sqrt(P+1)+sqrt(P)).^(2/3);
    t=(1+Q+1./Q)/6;
    %The max command prevents finite precision problems due to
    %subtraction within the square root.
    c=max(0,u.^2-1+2*t);
    c=sqrt(c);
    w=(c-u)/2;

    %The z coordinate of the closest point projected on the ellipsoid.
    %The max command deals with precision problems when the argument
    %is nearly zero. The problems arise due to the subtraction within
    %the square root.
    z=max(0,sqrt(t.^2+v)-u.*w-t/2-1/4);
    z=sign(z0).*sqrt(q).*(w+sqrt(z));
    Ne=a*sqrt(1+eps2.*z.^2/b^2);

    %The min and max terms deals with finite precision problems.
    val=min(z*(eps2+1)./Ne,1);
    val=max(val,-1);
    phi(curPoint)=asin(val);
    h(curPoint)=r0.*cos(phi(curPoint))+z0.*sin(phi(curPoint))-a^2./Ne;
end
    
end

function [lambda,phi,h]=FukishimaAlg(cartPoints,a,f)
%%FUKUSHIMAALG This function implements the algorithm of [1] with minor
%              modifications.
%
%If one lets the algorithm run for an arbitrary number of iterations, there
%will generally be underflows, since the ratio of S and C matter, but both
%terms can drift by a constant factor during the iterations. Thus, after
%each iterative step, the values are normalized so that C=1 (it takes one
%division). Convergence of Fukushima's method is assumed to occur after
%six iterations.
%
%REFERENCES:
%[1] Fukushima, T., "Transformation from Cartesian to geodetic coordinates
%    accelerated by Halley's method", J.Geodesy (2006) 79: 689-693.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPoints=size(cartPoints,2);

%The semi-minor axis of the reference ellipsoid.
b=a*(1-f);

%The square of the first numerical eccentricity. 
e2=2*f-f^2;

ec=sqrt(1-e2);

%Allocate space for the results.
phi=zeros(1,numPoints);
lambda=zeros(1,numPoints);
h=zeros(1,numPoints);
for curPoint=1:numPoints
    %Extract the coordinates
    x0=cartPoints(1,curPoint);
    y0=cartPoints(2,curPoint);
    z0=cartPoints(3,curPoint);
    
    lambda(curPoint)=atan2(y0,x0);

    p=sqrt(x0^2+y0^2);
    P=p/a;
    Z=(ec/a)*abs(z0);

    S=Z;
    C=ec*P;

    %Loop until convergence. Normally, only 6 iterations or fewer is
    %required. When some points are near the surface of the Earth and f is
    %close to 1, the required number of iterations can be higher.
    for curIter=1:500
        A=sqrt(S^2+C^2);
        B=1.5*e2*S*C^2*((P*S-Z*C)*A-e2*S*C);
        F=P*A^3-e2*C^3;
        D=Z*A^3+e2*S^3;

        SNew=D*F-B*S;
        CNew=F^2-B*C;
        
        SOld=S;
        COld=C;

        SNew=SNew/CNew;
        if(~isfinite(SNew))
            S=SNew;
            break;
        else
            S=SNew;
            C=1;
        end
        
        if(S==SOld&&C==COld)
           break;
        end
    end
    Cc=ec*C;

    %If the point is along the z-axis, then SNew and CNew will
    %both be zero, leading to a non-finite result.
    if(~isfinite(S))
        phi(curPoint)=sign(z0)*(pi/2);
        h(curPoint)=abs(z0)-b;
    else
        phi(curPoint)=sign(z0)*atan(S/Cc);
        h(curPoint)=(p*Cc+abs(z0)*S-b*A)/sqrt(Cc^2+S^2);
    end
end
    
end

function [lambda,phi,h]=OlsonAlg(cartPoints,a,f)
%%OLSONALG This function implements the algorithm of [1], removing a test
%          for a radius being too small as the function will still work at
%          smaller radii.
%
%REFERENCES:
%[1] D. K. Olson, "Converting Earth-centered, Earth-fixed coordinates to
%   geodetic coordinates," IEEE Transactions on Aerospace and Electronic
%   Systems, vol. 32, no. 1, pp. 473-476, Jan. 1996.
%
%October 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numPoints=size(cartPoints,2);

%The square of the eccentricity.
e2=2*f-f^2;

a1=a*e2;
a2=a1^2;
a3=a1*e2/2;
a4=(5/2)*a2;
a5=a1+a3;
a6=1-e2;

%Allocate space for the results.
phi=zeros(1,numPoints);
lambda=zeros(1,numPoints);
h=zeros(1,numPoints);
for curPoint=1:numPoints
    %Extract the coordinates
    x=cartPoints(1,curPoint);
    y=cartPoints(2,curPoint);
    z=cartPoints(3,curPoint);
    
    zp=abs(z);
    w2=x^2+y^2;
    w=sqrt(w2);
    z2=z^2;
    r2=w2+z2;
    %The algorithm will work with points close to the origin. Thus, there
    %is no need to have a test for r being too small as is the case in [1].
    r=sqrt(r2);

    lambda(curPoint)=atan2(y,x);
    s2=z2/r2;
    c2=w2/r2;
    u=a2/r;
    v=a3-a4/r;
    if(c2>0.3)
        s=(zp/r)*(1+c2*(a1+u+s2*v)/r);
        phi(curPoint)=asin(s);
        ss=s^2;
        c=sqrt(1-ss);
    else
        c=(w/r)*(1-s2*(a5-u-c2*v)/r);
        phi(curPoint)=acos(c);
        ss=1-c^2;
        s=sqrt(ss);
    end

    g=1-e2*ss;
    rg=a/sqrt(g);
    rf=a6*rg;
    u=w-rg*c;
    v=zp-rf*s;
    f=c*u+s*v;
    m=c*v-s*u;
    p=m/(rf/g+f);
    phi(curPoint)=phi(curPoint)+p;
    h(curPoint)=f+m*p/2;
    if(z<0)
        phi(curPoint)=-phi(curPoint);
    end
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
