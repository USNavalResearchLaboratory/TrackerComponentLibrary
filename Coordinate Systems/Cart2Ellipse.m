function points=Cart2Ellipse(cartPoints,a,f)
%%CART2ELLIPSE Convert Cartesian coordinates to ellipsoidal coordinates.
%
%INPUTS: cartPoints A matrix of the points in ECEF Cartesian coordinates
%                   that are to be transformed into ellipsoidal
%                   coordinates. Each column of cartPoints is of the
%                   format [x;y;z].
%           a       The semi-major axis of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84SemiMajorAxis is used.
%           f       The flattening factor of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84Flattening is used.
%
%OUTPUTS:   points  A matrix of the converted points. Each column of the
%                   matrix has the format [latitude;longitude;altitude],
%                   with latitude and longitude given in radians.
%
%The primary algorithm implemented in this function is [1], which is a
%modified version of the algorithm in [2].
%
%However, both techniques cited above will fail if the point in question
%is too close to the origin (deep within the Earth). In such an instance,
%Halley's method as described in [3] is used instead, with a minor
%modification. If one lets the algorithm run for an arbitrary number of
%iterations, there will generally be underflows, since the ratio of S and C
%matter, but both terms can drift by a constant factor during the
%iterations. Thus, after each iterative step, the values are normalized so
%that C=1 (it takes one division). Convergence of Fukushima's method is
%assumed to occur after six iterations.
%
%REFERENCES:
%[1] I. Sofair "Improved method for calculating exact geodetic latitude and
%    altitude revisited," Journal of Guidance, Control, and Dynamics, vol.
%    23, no. 2, p. 369, Mar. 2000.
%[2] I. Sofair, "Improved method for calculating exact geodetic latitude
%    and altitude," Journal of Guidance, Control, and Dynamics, vol. 20,
%    no. 4, pp. 824-826, Jul.-Aug. 1997.
%[3] Fukushima, T., "Transformation from Cartesian to geodetic coordinates
%    accelerated by Halley's method", J.Geodesy (2006) 79: 689-693.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    f=Constants.WGS84Flattening;
end

if(nargin<2)
    a=Constants.WGS84SemiMajorAxis;
end

numPoints=size(cartPoints,2);

%The semi-minor axis of the reference ellipsoid.
b=a*(1-f);

%The square of the first numerical eccentricity. 
e2=2*f-f^2;
%The square of the second numerical eccentricity.
eps2=a^2/b^2-1;

%This value is used if Fukushima's method is chosen.
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
    
    r0=sqrt(x0.^2+y0.^2);
    p=abs(z0)/eps2;
    s=r0.^2/(e2*eps2);
    q=p.^2-b.^2+s;
    
    lambda(curPoint)=atan2(y0,x0);
    
    if(q>=0)%Use Sofair's algorithm
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
    else%Use Fukushima's algorithm.
        p=sqrt(x0^2+y0^2);
        P=p/a;
        Z=(ec/a)*abs(z0);
        
        S=Z;
        C=ec*P;
        
        %Loop until convergence. Assume convergence in 6 iterations.
        for curIter=1:6
            A=sqrt(S^2+C^2);
            B=1.5*e2*S*C^2*((P*S-Z*C)*A-e2*S*C);
            F=P*A^3-e2*C^3;
            D=Z*A^3+e2*S^3;

            SNew=D*F-B*S;
            CNew=F^2-B*C;

            SNew=SNew/CNew;
            if(~isfinite(SNew))
                S=SNew;
                break;
            else
                S=SNew;
                C=1;
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

points=[phi;lambda;h];

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
