function latLonHOsc=ellips2SpecOscCoords(latLonHPts,rOsc,spherCent,aEllips,fEllips)
%%ELLIPS2SPECOSCCOORDS Given points in [latitude;longitude;height]
%           ellipsoidal coordinates, convert them into the
%           [latitude;longitude;height] system of a specific sphere
%           (usually an osculating sphere is the most useful) defined with
%           radius rOsc and center spherCent with respect to the global
%           coordinate system.
%
%INPUTS: latLonHPts The 3XN set of [latitude;longitude;height] points to
%               convert from ellipsoidal to coordinates on the specified
%               sphere. Latitude and longitude are given in radians.
%               If a 2XN matrix is passed, it is assumed that the omitted
%               height coordinate is zero and the height is omitted from
%               latLonHOsc.
%          rOsc The scalar radius of the osculating sphere.
%    sphereCent The 3X1 location of the osculating sphere in the global
%               coordinate system where the origin is the center of the
%               reference ellipsoid.
%       aEllips The semi-major axis of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84SemiMajorAxis is used.
%       fEllips The flattening factor of the reference ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84Flattening is used.
%
%OUTPUTS: latLonHOsc A 3XN set of points in [latitude;longitude;height] in
%                    the osculating sphere coordinate system. If latLonHPts
%                    was 2XN, then the height component of this matrix is
%                    also omitted.
%
%The use of osculating spheres is discussed in [1].
%
%REFERENCES:
%[1] P. Williams and D. Last, "On Loran-C time-difference to co-ordinate
%    converters," in Proceedings of the 32nd Annual Convention & Technical
%    Symposium of the International Loran Association, Boulder, CO, 3-7
%    Nov. 2003.
%
%December 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<5||isempty(fEllips))
        fEllips=Constants.WGS84Flattening;
    end

    if(nargin<4||isempty(aEllips))
        aEllips=Constants.WGS84SemiMajorAxis;
    end

    ptsCartEllips=ellips2Cart(latLonHPts,aEllips,fEllips);
    ptsCartSpher=bsxfun(@minus,ptsCartEllips,spherCent);
    latLonHOsc=Cart2Ellipse(ptsCartSpher,[],rOsc,0);
    
    if(size(latLonHPts,1)==2)
        %If the height should be omitted.
        latLonHOsc=latLonHOsc(1:2,:);
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
