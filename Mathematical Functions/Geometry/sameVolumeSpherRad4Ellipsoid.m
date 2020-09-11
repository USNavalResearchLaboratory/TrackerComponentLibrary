function rVol=sameVolumeSpherRad4Ellipsoid(a,f)
%%SAMEVOLUMESPHER4ELLIPSOID Given an ellipsoid of rotation, determine the
%               radius of a sphere having the same surface volume.
%
%INPUTS: a The semi-major axis length of the ellipsoid. If this argument is
%          omitted or an empty matrix is passed, the value in
%          Constants.WGS84SemiMajorAxis is used.
%        f The flattening factor of the ellipsoid. If this argument is
%          omitted or an empty matrix is passed, the value in
%          Constants.WGS84Flattening is used.
%
%OUTPUTS: rVol The radius of a sphere having the same volume as the
%              specified ellipsoid.
%
%An expression for ther same-volume sphere is given in [1].
%
%REFERENCES:
%[1] H. Mortiz, "Geodetic Reference System 1980," Bulletin Géodésique,
%    vol. 54, no. 3, pp. 395-405, Sep. 1980. Given with corrections at
%    https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<1||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

%The semi-minor axis of the ellipsoid.
b=a*(1-f);

rVol=nthroot(a^2*b,3);

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
