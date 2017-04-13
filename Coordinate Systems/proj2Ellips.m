function zGround=proj2Ellips(z,a,f)
%%PROJ2ELLIPS   Given a point in three-dimensional space, find the point on
%               a reference ellipsoid that is closest. This approximately
%               finds the point on the ground (ignoring terrain) that is
%               closest to a given arbitrary point.
%
%INPUTS:    z       A point in [x;y;z] Cartesian coordinates.
%           a       The semi-major axis of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84SemiMajorAxis is used.
%           f       The flattening factor of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84Flattening is used.
%OUTPUTS:  zGround  The point on the surface of the reference ellipsoid
%                   that is closest to the ground.
%
%This just transforms the given point into latitude, longitude and
%altitude, sets the altitude to zero and converts back to Cartesian
%coordinates. The function used for the conversion from Cartesian to
%ellipsoidal coordinates assumes that the point in question is not too
%close to the center of the Earth.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<2||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

zEllipse=Cart2Ellipse(z,[],a,f);
zEllipse(3)=0;
zGround=ellips2Cart(zEllipse,a,f);
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
