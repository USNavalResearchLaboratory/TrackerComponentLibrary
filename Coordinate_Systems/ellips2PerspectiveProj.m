function xyh=ellips2PerspectiveProj(plhPoints,plhCamera,a,f)
%%ELLIPS2PERSPECTIVEPROJ Convert points given in ellipsoidal coordinates
%       over the Earth into a perspective projection with respect to a
%       particular camera location. The camera is facing straight down at
%       the Earth.
%
%INPUTS: plhPoints A 3XN set of [latitude;longitude;height] points to
%               convert with latitude and longitude given in radians.
%     plhCamera The 3X1 [latitude;longitudeheight] reference point for the
%               camera. h must be finite.
%             a The semi-major axis of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84Flattening is used.
%
%OUTPUTS: xyh The 3XN set of [x;y;h] points of the perspective projection.
%             The h components are the same as the heights on the input.
%
%This implements the general perspective projection for an ellipsoid that
%is given in Chapter 23 of [1].
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<3||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

e2=f*(2-f);%Squared eccentricity.

if(~isfinite(plhCamera(3)))
    P=Inf;
else
    lCamera=ellips2Cart(plhCamera,a,f);
    P=norm(lCamera)/a;
end

%Camera coordinates:
phi1=plhCamera(1);
lambda0=plhCamera(2);
H=plhCamera(3);

%Coordinates of the points.
phi=plhPoints(1,:);
lambda=plhPoints(2,:);
h=plhPoints(3,:);

sinPhi=sin(phi);
sinPhi1=sin(phi1);
cosPhi=cos(phi);
cosPhi1=cos(phi1);
sinLambdaDiff=sin(lambda-lambda0);
cosLambdaDiff=cos(lambda-lambda0);

%Equation 4-20
N=a./sqrt(1-e2*sinPhi.^2);
%Equationa/(1-e2 8-23
N1=a./sqrt(1-e2*sinPhi1.^2);
%Equation 23-15
C=((N+h)./a).*cosPhi;

%Equation 23-16
S=((N.*(1-e2)+h)./a).*sinPhi;

%Equation 23-17
phig=phi1-asin(N1.*e2.*sinPhi1.*cosPhi1./(P.*a));
%Equation 23-19
K=H./(P.*cos(phi1-phig)-S.*sinPhi1-C.*cosPhi1.*cosLambdaDiff);
%Equation 23-19a and Equation 23-20
xyh=[K.*C.*sinLambdaDiff;
    K.*(P.*sin(phi1-phig)+S.*cosPhi1-C.*sinPhi1.*cosLambdaDiff);
    h];
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


