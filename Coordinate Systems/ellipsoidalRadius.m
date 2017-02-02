function r=ellipsoidalRadius(el,elType,a,f)
%%ELLIPSOIDALRADIUS Find the distance from the center of an axis-aligned
%           oblate spheroid, as is used to approximate the shape of the
%           Earth, to a point on the surface of the spheroid as a function
%           of the (spherical) elevation angle above the equatorial plane
%           or as a function of the ellipsoidal latitude.
%
%INPUTS: el A point or a matrix of the elevation points for which one wants
%           to find radii. The points are in radians. The next parameter
%           sets the type of points.
%    elType An optional parameter specifying the type of elevation el is.
%           Possible values are:
%           0 (The default if omitted or an empty matrix is passed)
%             Geocentric elevations (as in spherical cooridinates).
%           1 Geodetic (ellipsoidal) latitudes.
%         a The semi-major axis of the reference ellipsoid. If this
%           argument is omitted, the value in Constants.WGS84SemiMajorAxis
%           is used.
%         f The flattening factor of the reference ellipsoid. If this
%           argument is omitted, the value in Constants.WGS84Flattening is
%           used.
%
%OUTPUTS: r A matrix of the radii corresponding to the values in el.
%
%The formula for the radius as a function of spherical elevation is taken
%from Section 3.2 of [1].
%
%REFERENCES:
%[1] C. Hirt and M. Rexer, "Earth2014: 1 arc-min shape, topography, bedrock
%    and ice-sheet models - available as gridded data and degree - 10,800 
%    spherical harmonics," International Journal of Applied Earth
%    Observation and Geoinformation, vol. 39, pp. 103-112, Jul. 2015.
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(elType))
   elType=0; 
end

if(nargin<3||isempty(a))
   a=Constants.WGS84SemiMajorAxis;
end

if(nargin<4||isempty(f))
   f=Constants.WGS84Flattening;
end

b=a*(1-f);%Semi-minor axis
e2=(a^2-b^2)/a^2;%eccentricity squared

if(elType==0)
    %We need the square of the sine of the ellipsoidal latitude. However,
    %we only have the spherical elevation. Here, we just substituted the
    %formula from spherLat2EllipsLat into the square of the sine and
    %simplified to get the following formula in terms of the squared
    %tangent.
    spherLat=wrapRange(el,-pi/2,pi/2,true);
    sinEl2=1-1./(1+(1/(1-e2))^2*tan(spherLat).^2);
else
    sinEl2=sin(el).^2;
end

r=a*sqrt((1-e2*(2-e2)*sinEl2)./(1-e2*sinEl2));

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
