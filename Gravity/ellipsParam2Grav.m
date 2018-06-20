function [U,g]=ellipsParam2Grav(cartPoint,omega,a,f,GM)
%%ELLIPSPARAM2GRAV Obtain the gravity potential and acceleration on
%                  the surface of the reference ellipsoid (this potential
%                  defines the geoid) or at an arbitrary point, and/ or get
%                  the acceleration due to gravity (INCLUDING centrifugal
%                  forces due to the Earth's rotation) that are
%                  associated with the defining ellipsoid parameters. The
%                  gravitational potential is the potential energy per
%                  unit mass.
%
%INPUTS: cartPoint The 3XN matrix of N points at which the potential and
%                  gravitational acceleration are desired given in ECEF
%                  [x;y;z] coordinates. If this parameter is omitted, or if
%                  an empty matrix is passed, then only the gravitational
%                  potential on the surface of the ellipsoid will be
%                  returned. That potential does not depend on the
%                  location.
%            omega The average rotation rate of the planet. If this
%                  argument is omitted, the value in
%                  Constants.WGS84EarthRotationRate is used, having units
%                  of radians per second.
%                a The semi-major axis of the reference ellipsoid. If this
%                  argument is omitted, the value in
%                  Constants.WGS84SemiMajorAxis is used, having units of
%                  meters.
%                f The flattening factor of the reference ellipsoid. If
%                  this argument is omitted, the value in
%                  Constants.WGS84Flattening is used (and is unitless).
%               GM The universal gravitational constant times the mass of
%                  the planet. If this parameter is omitted, then the value
%                  in Constants.WGS84GMWithAtmosphere is used, having units
%                  of m^3/s^2.
%
%OUTPUTS: U The NX1 vector of gravity potentials at the desired points or
%           on the surface of the reference ellipsoid if no point is given.
%           The gravity potential has units of m^2/s^2 and includes the
%           effects of the rotation of the Earth. 
%         g The 3XN matrix of acceleration vectors due to gravity at the
%           desired points in Cartesian coordinates, including the
%           acceleration due to the Earth's rotation. This is the gradient
%           of U and points toward the Earth. This parameter cannot be
%           returned if cartPoint is empty or omitted.
%
%This function relates the defining parameters for a reference ellipsoid,
%such as that in the DoD's WGS-84 standard, to a gravitational potential on
%the surface of the ellipsoid.
%
%The formula for the potential on the surface of the reference ellipsoid
%is Equation 2-123 in Chapter 2-7 of [1]. The formula for the gravity
%potential at an arbitrary point is Equation 2-126. The formula for the
%acceleration due to gravity is Equation 2-132 in Chapter 2.8, which is
%given in terms of ellipsoidal parameters that are then converted to
%Cartesian coordinates.
%
%REFERENCES:
%[1] B. Hofmann-Wellenhof and H. Moritz, Physical Geodesy, 2nd ed. 
%    SpringerWienNewYork, 2006.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5)
    GM=Constants.WGS84GMWithAtmosphere;
end

if(nargin<4)
    f=Constants.WGS84Flattening;
end

if(nargin<3)
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<2)
    omega=Constants.WGS84EarthRotationRate;
end

%The first numerical eccentricity.
b=a*(1-f);
E=sqrt(a^2-b^2);

if(nargin<1||isempty(cartPoint))
    if(nargout==2)
        error('The gravitational acceleration is requested, but no point was provided')
    else
    %Just the geoid potential is desired. This is Equation 2-123.
        U=GM/E*atan(E/b)+(1/3)*omega^2*a^2;
        return;
    end
end

numPoints=size(cartPoint,2);

%Convert the point to ellipsoidal harmonic coordinates.
harmonCoord=Cart2EllipsHarmon(cartPoint,E);
beta=harmonCoord(1,:)';
u=harmonCoord(3,:)';

w=sqrt((u.^2+E^2*sin(beta).^2)./(u.^2+E^2));
qp=3*(1+u.^2/E^2).*(1-u/E.*atan(E./u))-1;
q0=0.5*((1+3*b^2/E^2)*atan(E/b)-3*b/E);
q =0.5*((1+3*u.^2/E^2).*atan(E./u)-3*u/E);

%First, find the gravity potential. This is Equation 2-126.
U=(GM/E)*atan(E./u)+(1/2)*omega^2*a^2*(q./q0).*(sin(beta).^2-1/3)+(1/2)*omega^2*(u.^2+E^2).*cos(beta).^2;

%Next, find the components of the acceleration vector. These are from
%Equation 2-132. If E==0, then the q values will be NaNs and we will
%substitute asymptotic values.
qRat=(q./q0);
qpRat=(qp./q0);
gammaLambda=zeros(numPoints,1);
if(any(~isfinite(qRat))||any(~isfinite(qpRat)))
    qRat=b^3./u.^3;
    gammau=-(GM./(u.^2)+(omega^2*a^2*3*b^3./u.^4).*(0.5*sin(beta).^2-1/6)-omega^2*u.*cos(beta).^2)./w;
else
    gammau=-(GM./(u.^2+E^2)+(omega^2*a^2*E./(u.^2+E^2)).*qpRat.*(0.5*sin(beta).^2-1/6)-omega^2*u.*cos(beta).^2)./w;
end
gammaBeta=-((-omega^2*a^2./sqrt(u.^2+E^2).*qRat+omega^2*sqrt(u.^2+E^2)).*sin(beta).*cos(beta))./w;

%Get the direction vectors for the given ellipsoidal harmonic coordinates.
g=zeros(3,numPoints);
for curPoint=1:numPoints
    u=getEllipsHarmAxes(harmonCoord(:,curPoint),E);

    %Find the gravitational acceleration vector in Cartesian space. The
    %negative sign is because the acceleration due to gravity
    g(:,curPoint)=u(:,1)*gammaBeta(curPoint)+u(:,2)*gammaLambda(curPoint)+u(:,3)*gammau(curPoint);
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
