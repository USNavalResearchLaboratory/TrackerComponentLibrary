function [xyZ,sigma,gammaVal]=ellips2PolarStereo(latLon,lambda0,k0,xyPole,a,f)
%%ELLIPS2POLARSTEREO Convert a point on the reference ellipsoid in
%            ellipsoidal (geodetic) coordinates gives as a latitude and a
%            longitude into a point on the polar stereographic projection.
%            Stereographic projections tend to be used near the poles. The
%            Uniform Polar Stereographic (UPS) coordinate system, which is
%            a specific realization of the polar stereographic coordinate
%            system, is typically only used at latitudes >=84 degrees and
%            those < -80 degrees.
%
%INPUTS: latLon A 2XN set of N [latitude;longitude] pairs given in radians
%               to be converted into polar stereographic coordinates.
%       lambda0 The location in raidns of the central meridian. If this
%               parameter is omitted or an empty matrix is passed, the
%               default of 0 is used, which is the value used in the UPS
%               system.
%            k0 The unitless central scale factor. If this parameter is
%               omitted or an empty matrix is passed, then the default of
%               0.994 is used, which is the value used in the UPS.
%        xyPole A 2X1 vector holding the false Easting and False
%               Northing of the origin. This is typically given in meters.
%               If this parameter is omitted or an empty matrix is passed,
%               then the default of [2000000;2000000] is used, which is the
%               value in meters used for the UPS system.
%             a The semi-major axis of the reference ellipsoid. If this
%               argument is omitted, the value in
%               Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted, the value in Constants.WGS84Flattening
%               is used.
%
%OUTPUTS: xyZ A 3XN set of points in polar stereographic coordinates
%             consisting of [Easting;Northing;Z], where Z=sign(latitude),
%             so 1 for the northern hemisphere and -1 for the southern
%             hemisphere and 0 on the equator. Easting and Northing will be
%             in meters if a and xyPole are in meters.
%       sigma The NX1 set of point scales. This indicates how the map
%             projections enlarges or reduces small distances in this type
%             of map projection.
%    gammaVal The NX1 set of convergence of meridians in radians. This
%             gives the angles of intersection between the meridian and
%             vertical lines on the map projection that are given by
%             x=constant value. This value
%
%If lambda0, k0, xyPole, a and f are all omitted, then the conversion is
%for the UPS on the WGS-84 reference ellipsoid.
%
%Chapter 21 of [1] discusses the polar stereographic projection. The
%implementation is taken from Sections 8 and 9 of [1]. Section 10 of 2
%provides the values for the UPS.
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%[2] Office of Geomatics, "National geospatial-intelligence agency
%    standardization document: Implementation practice: The universal
%    grids and the transverse mercator and polar stereographic map
%    projections," National Geospatial-Intelligence Agency, Tech. Rep.
%    NGA.SIG.0012 2.0.0 UTMUPS, 25 Mar. 2014. [Online]. Available:
%    http://earth-info.nga.mil/GandG/publications/NGA_SIG_0012_2_0_0_UTMUPS/NGA.SIG.0012_2.0.0_UTMUPS.pdf
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(lambda0))
   lambda0=0;%Central meridian location for UPS.
end

if(nargin<3||isempty(k0))
   k0=0.994;%Central scale factor for UPS
end

if(nargin<4||isempty(xyPole))
   xyPole=[2000000;2000000];%Easting and Northing of the pole for UPS.
end

if(nargin<5||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<6||isempty(f))
    f=Constants.WGS84Flattening;
end

%The first numerical eccentricity of the ellipsoid.
e=sqrt(2*f-f^2);

N=size(latLon,2);

xyZ=zeros(3,N);
sigma=zeros(N,1);
gammaVal=zeros(N,1);
for curPoint=1:N
    phi=latLon(1,curPoint);
    lambda=latLon(2,curPoint);
    lambda=lambda-lambda0;

    Z=sign(phi);
    if(Z==-1)
        phi=-phi;
    end

    sinPhi=sin(phi);
    sinLambda=sin(lambda);
    cosLambda=cos(lambda);

    [sinChi,cosChi]=ellipsLat2SinCosConformLat(phi,f);

    k90=sqrt(1-e^2)*exp(e*atanh(e));
    denom=k90*(1+sinChi);
    x=2*k0*a*sinLambda*cosChi/denom;
    y=-2*k0*a*cosLambda*cosChi/denom;
    sigma(curPoint)=k0*2*sqrt(1-e^2*sinPhi^2)*exp(e*atanh(e*sinPhi))/(k90*(1+sinPhi));
    gammaVal(curPoint)=lambda;

    if(Z==-1)
        y=-y;
        gammaVal(curPoint)=-gammaVal(curPoint);
    end

    xyZ(:,curPoint)=[x;y;Z]+[xyPole;0];
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
