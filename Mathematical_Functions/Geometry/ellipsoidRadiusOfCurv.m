function RAlpha=ellipsoidRadiusOfCurv(phi,alpha,a,f)
%%ELLIPSOIDRADIUSOFCURV Find the radius of curvature at a particular point
%           in a particular direction on the surface of an axis-aligned
%           ellipsoid. That is, we are finding the radius of curvature in
%           the normal section azimuth alpha. The semi-major axis falls in
%           the x-y plane and the semi-minor axis is in the z direction
%           (North).
%
%INPUTS: phi The latitude in radians where the radius of curvature is
%            desired.
%      alpha The azimuth direction in radians Easth of North. In other
%            words, one wishes to get the local radius of curvature for
%            someone traveling in this direction. alpha=0 corresponds to
%            the radius of curvature in the meridian (meridians pass
%            through the poles) and alpha=+/-pi/2 corresponds to the radius
%            of curvature in the prime vertical. That is considering the
%            curvature of the curve resulting from the cut of the ellipsoid
%            that is normal to the meridian at a point as well as normal to
%            the surface of the ellipsoid at that point.
%          a The semi-major axis length of the ellipsoid. If this argument
%            is omitted or an empty matrix is passed, the value in
%            Constants.WGS84SemiMajorAxis is used.
%          f The flattening factor of the ellipsoid. If this argument is
%            omitted or an empty matrix is passed, the value in
%            Constants.WGS84Flattening is used.
%
%OUTPUTS: RAlpha The desired scalar radius of curvature.
%
%Equation 3.104 in [1] implements the radius of curvature in the normal
%section azimuth alpha.
%
%EXAMPLE:
%Here, we show that the values for 0 and pi/2 radius respectively equal the
%radius of curvature in the meridian and the radius of cuavture in the
%prime vertical.
% phi=deg2rad(26);
% f=Constants.WGS84Flattening;
% a=Constants.WGS84SemiMajorAxis;
% e2=f*(2-f);%The squared eccentricity.
% %Equation 3.17 in [1]. The radius of curvature in the meridian.
% M=a*(1-e2)/sqrt(1-e2*sin(phi)^2)^3;
% alpha=0;
% MAlpha=ellipsoidRadiusOfCurv(phi,alpha,a,f);
% %The two values are equal withing finite precision limits.
% norm(M-MAlpha)
% %Equation 3.99 in [1]. The radius of curvature in the prime vertical.
% N=a/sqrt(1-e2*sin(phi)^2);
% alpha=pi/2;
% NAlpha=ellipsoidRadiusOfCurv(phi,alpha,a,f);
% %The two values are equal within finite precision limits.
% norm(NAlpha-N)
%
%REFERENCES
%[1] R. H. Rapp, "Geometric geodesy, part I," Ohio State University
%    Department of Geodetic Science and Surveying, Tech. Rep., Apr.
%    1991. [Online]. Available: http://hdl.handle.net/1811/24333
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

if(nargin<2||isempty(alpha))
    alpha=0;
end

%The squared eccentricity of the ellipsoid.
e2=f*(2-f);
ePrime2=e2/(1-e2);

%The radius of curvature in the prime vertical. This considers the curve
%resulting from the cut of the ellipsoid that is normal to the meridian at
%a point as well as normal to the surface of the ellipsoid at that point.
%Equation 3.99 in [1].
N=a/sqrt(1-e2*sin(phi)^2);

%Equation 3.104.
RAlpha=N/(1+ePrime2*cos(alpha)^2*cos(phi)^2);

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
