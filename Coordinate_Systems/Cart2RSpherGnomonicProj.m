function [rxGyG,cosC]=Cart2RSpherGnomonicProj(xyzPts,latLonRef)
%%CART2RSPHERGNOMONICPROJ Convert Cartesian points (given with respect to a
%          local reference sphere) to a Gnomonic projection defined at a
%          specific spherical latitude and longitude. The assumed spherical
%          radius of the conversion is returned as the first coordinate of
%          each transformed point, which differs from spher2GnomonicProj
%          where the radius is specified and points are only addressed on
%          the reference sphere surface by latitude and longitude. This is
%          the projection from the Cartesian points on a sphere to the
%          plane such that a line connecting the point on the sphere and
%          the point on the plane is normal to the sphere.
%
%INPUTS: xyzPts A 3XN set of points on the reference sphere in [x;y;z]
%               Cartesian coordinates that should be converted to
%               gnomonically projected coordinates in the tangent plane.
%     latLonRef The 2X1 [latitude;longitude] reference point for the
%               gnomonic projection.
%
%OUTPUTS: rxGyG A 3XN set of the points gnomonically projected into the
%            local tangent plane at latLonRef. The coordinates are
%            [r;xG;yG], whereby xG and yG are the gnomonic coordinates and
%            r is the radius of the spherical approximation. The 3D
%            location of the ith point in the plane is
%            rxGyG(2,i)*uEast+rxGyG(3,i)*uNorth, where uEast and uNorth are
%            the first two vectors in the output of
%            getENUAxes(latLonRef,false,rxGyG(1,i),0).
%       cosC This is the 1XN seet of values cosine of the distance between
%            points in xyzPts and latLonRef when taken on the surface of
%            the unit sphere (distance is radians). If value of this are
%            negative, then the point in question is actually on the other
%            side of the Earth from the reference point, so the inverse
%            projection gnomonicProj2Sphere will put it on the wrong side
%            of the Earth.
%
%The radius of the reference sphere association with each point is taken to
%be the magnitude of the point. The conversion comes from substituting
%Cartesian relations into Chapter 22 of [1].
%
%EXAMPLE:
%In this instance, spher2GnomonicProj is used to get the gnomonic
%coordinates of a point that has been specified in terms of latitude and
%longitude. That point is then converted to Cartesian coordinates and this
%function is used to get the same Gnomonic coordinates. The relative error
%between the conversions is on the order of finite precision limits, which
%indicates that the functions are consistent with each other.
% latLonRef=deg2rad([20.756113;-156.010933]);
% latLonPt=deg2rad([19.456233;-154.823617]);
% r=osculatingSpher4LatLon(latLonRef);
% xyGnomon=spher2GnomonicProj(latLonPt,latLonRef,r);
% xyzCart=ellips2Cart(latLonPt,r,0);
% rxyGnomon=Cart2RSpherGnomonicProj(xyzCart,latLonRef);
% %The relative error between the direct spherical conversion and the
% %Cartesian conversion.
% RelErr=max(abs((xyGnomon-rxyGnomon(2:3,:))./xyGnomon))
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The sine and cosine of the reference latitude.
sinPhi1=sin(latLonRef(1));
cosPhi1=cos(latLonRef(1));
%The sine and cosine of the reference longitude.
sinLambda0=sin(latLonRef(2));
cosLambda0=cos(latLonRef(2));

x=xyzPts(1,:);
y=xyzPts(2,:);
z=xyzPts(3,:);
r=sqrt(sum(xyzPts.*xyzPts,1));

rxGyG=[r;
       (-r.*(x.*sinLambda0-y.*cosLambda0))./(x.*cosPhi1.*cosLambda0+y.*cosPhi1.*sinLambda0+z.*sinPhi1);
       (-r.*(x.*sinPhi1.*cosLambda0+y.*sinPhi1.*sinLambda0-z.*cosPhi1))./(x.*cosPhi1.*cosLambda0+y.*cosPhi1.*sinLambda0+z.*sinPhi1)];

if(nargout>1)
    cosC=(x.*cosPhi1.*cosLambda0+y.*cosPhi1.*sinLambda0+z.*sinPhi1)./r;
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
