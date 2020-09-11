function uDot=uDotEllipsoid(u,x,a,f)
%%UDOTELLIPSOID Return the derivative of the basis vectors for the local
%               coordinate system of a moving observer when using an
%               ellipsoidal approximation of the curved Earth. Such an
%               evolving coordinate system has been termed "wander
%               coordinate" or "naturally evolving coordinates". Targets
%               moving locally in non-maneuvering, level flight with
%               a coordinate system evolving as per this coordinate system
%               will follow geodesic curves, not rhumb lines.
%
%INPUTS: u The 3X3 orthonormal basis vectors for the local coordinate
%          system. u(:,3) is the local vertical.
%        x An nX1 vector whose first three elements are Cartesian position
%          in the global ECEF coordinate system and whose next three
%          elements are velocity in the local (flat-Earth) coordinate
%          system. Other elements of x do not matter.
%        a The semi-major axis of the reference ellipsoid. If this argument
%          is omitted, the value in Constants.WGS84SemiMajorAxis is used.
%        f The flattening factor of the reference ellipsoid. If this
%          argument is omitted, the value in Constants.WGS84Flattening is
%          used.
%
%OUTPUTS: uDot The derivative of the basis vector u with respect to time.
%
%The solution is detailed in [1], and has been simplified for an
%ellipsoidal Earth so that numerical differentiation is not necessary.
%
%REFERENCES:
%[1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%    Earth's curvature," Journal of Advances in Information Fusion, vol.
%    10, no. 1, Jun. 2015.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<4||isempty(f))
        f=Constants.WGS84Flattening;
    end

    if(nargin<3||isempty(a))
        a=Constants.WGS84SemiMajorAxis;
    end

    %The first numerical eccentricity of the ellipsoid.
    e=sqrt(2*f-f^2);
    
    %The current position.
    r=x(1:3,1);
    %The local velocity.
    rDotLocal=x(4:6,1);
    
    %Get the velocity in global coordinates.
    rDot=getGlobalVectors(rDotLocal,u);
    
    %Convert to ellipsoidal coordinates.
    rEllipse=Cart2Ellipse(r,[],a,f);
    
    %The ellipsoidal latitude
    phi=rEllipse(1);
    
    %The height above the reference ellipsoid.
    h=rEllipse(3);
    
    %Obtain the local East-North-Up axes.
    ENUAxes=getENUAxes(rEllipse,false,a,f);
    EAxis=ENUAxes(:,1);
    NAxis=ENUAxes(:,2);
    
    %The local North and East components of the velocity vector.
    vE=dot(EAxis,rDot(1:3,1));
    vN=dot(NAxis,rDot(1:3,1));
 
    %The normal radius of curvature.
    Rp=a/sqrt(1-e^2*sin(phi)^2);
    
    %The radius of curvature of the meridian line.
    Rm=Rp*(1-e^2)/(1-e^2*sin(phi)^2);

    %The rotation rate North.
    omegaN=vE/(Rp+h);
    
    %The rotation rate East.
    omegaE=-vN/(Rm+h);
    
    %The 3D rotation vector.
    Omega=omegaN*NAxis+omegaE*EAxis;
    
    %Compute the derivatives of the coordinate system basis vectors with
    %respect to time.
    uDot(:,1)=cross(Omega,u(:,1));
    uDot(:,2)=cross(Omega,u(:,2));
    uDot(:,3)=cross(Omega,u(:,3));
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
