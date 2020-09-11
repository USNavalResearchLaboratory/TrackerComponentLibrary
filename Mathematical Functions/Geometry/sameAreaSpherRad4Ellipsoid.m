function rAuthalic=sameAreaSpherRad4Ellipsoid(a,f)
%%SAMEAREASPHERRAD4ELLIPSOID Given an ellipsoid of rotation, determine the
%               radius of a sphere having the same surface area. This is
%               known as the authalic radius. 
%
%INPUTS: a The semi-major axis length of the ellipsoid. If this argument is
%          omitted or an empty matrix is passed, the value in
%          Constants.WGS84SemiMajorAxis is used.
%        f The flattening factor of the ellipsoid. If this argument is
%          omitted or an empty matrix is passed, the value in
%          Constants.WGS84Flattening is used.
%
%OUTPUTS: rAuthalic The authalic radius for the specified 3D ellipsoid.
%
%As noted on page 4 of [1], authalic is from the Greek "autos" (same) and
%"ailos" (area). An expression for the authalic radius is on page 16 of
%[1],
%
%REFERENCES:
%[1] J. P. Snyder, Map Projections- A Working Manual, U. S. Geological
%    Survery Professional Paper 1935, United States Government Printing
%    Office: Washington, 1987.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<1||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

%The eccentricity of the ellipsoid.
e=sqrt(f*(2-f));

%Equation 3-13 in [1].
rAuthalic=(a/2)*sqrt(2+(e-1/e)*log(2/(1+e)-1));
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
