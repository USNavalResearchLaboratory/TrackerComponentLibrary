function pazPts=Cart2PolarAzEquidistProj(zCart,latLonRef,a,f,algSel)
%%CART2POLARAZEQUIDISTPROJ Convert 3D Cartesian points into an azimuthal
%       equidistant projection in polar form (ground range and azimuth
%       instead of Easting and Northing). These coordinates can be useful
%       for parameterizing tables of traced rays.
%
%INPUTS: xCart A 3XN set of Cartesian locations in an Earth centered Earth
%              fixed coordiante system.
%    latLonRef A 2X1 [latitude;longitude] reference point in radians about
%              which the projection is taken.
%            a The semi-major axis of the reference ellipsoid (in meters).
%              If this argument is omitted or an empty matrix is passed,
%              the value in Constants.WGS84SemiMajorAxis is used.
%            f The flattening factor of the reference ellipsoid. If this
%              argument is omitted or an empty matrix is passed, the value
%              in Constants.WGS84Flattening is used.
%       algSel An optional parameter specifying how the conversion is
%              performed. Possible values are:
%              0 (The default if omitted or an empty matrix is passed) Just
%                call the Cart2Ellipse function and then the
%                ellips2PolarAzEquidistProj function.
%              1 Use a simple conversion that is only valid for f=0 and
%                that is based on Equation 25-1, 5-3, and 5-4b in [1] and
%                substituting in the conversion from Cartsian to spherical
%                coordinates with an additional term for a height above the
%                sphere.
%
%OUTPUTS: pazPts A 3XN set of [ground distance; heading;height] points
%                where the heading is specified in radians East of North.
%
%EXAMPLE:
%This example demonstrates that the conversion is consistent with the
%polarAzEquidistProj2Cart when using algSel=0 using an osculating sphere 
%approximation for the Earth. A Cartesian point is converted into
%polar azimuthal coordinates and is then converted back. The relative error
%btwen the reverse conversion and the original Cartesian point is on the
%order of finite precision limitiations.
% latLonRef=deg2rad([20.906029;-157.078918]);
% a=osculatingSpher4LatLon(latLonRef);
% f=0;
% algSel=1;
% llhPt=[[deg2rad([20.631756;-155.551818]);0],[deg2rad([21.295516;-158.002386]);10e3]];
% xCart=ellips2Cart(llhPt,a,f);
% azEqCoords=Cart2PolarAzEquidistProj(xCart,latLonRef,a,f,algSel);
% xCartBack=polarAzEquidistProj2Cart(azEqCoords,latLonRef,a,f);
% RelErr=max(max(abs((xCartBack-xCart)./xCart)))
%
%REFERENCES:
%[1] J. P. Snyder, Map Projection - A Working Manual. Washington, D.C.:
%    U.S. Geological Survey, 1987, Professional Paper 1395.
%
%January 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(algSel))
    algSel=0; 
end

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

switch(algSel)
    case 0
        latLonPts=Cart2Ellipse(zCart,[],a,f);
        pazPts=ellips2PolarAzEquidistProj(latLonPts,latLonRef,a,f);
    case 1
        if(f~=0)
            error('Algorithm 1 only works with f=0')
        end
        
        x=zCart(1,:);
        y=zCart(2,:);
        z=zCart(3,:);

        rE=a;
        phi0=latLonRef(1);
        lambda0=latLonRef(2);

        sinPhi0=sin(phi0);
        cosPhi0=cos(phi0);
        sinLambda0=sin(lambda0);
        cosLambda0=cos(lambda0);

        r=sqrt(x.^2+y.^2+z.^2);
        pazPts=[rE.*acos((x.*cosLambda0.*cosPhi0+y.*cosPhi0.*sinLambda0+z.*sinPhi0)./r);
                 atan2(y.*cosLambda0-x.*sinLambda0,z.*cosPhi0-(x.*cosLambda0+y.*sinLambda0).*sinPhi0);
                 r-rE];
    otherwise
        error('Unknown value of algSel specified.')
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
