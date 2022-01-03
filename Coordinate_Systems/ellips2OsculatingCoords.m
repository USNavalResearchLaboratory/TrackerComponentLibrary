function [latLonHOsc,r,spherCent]=ellips2OsculatingCoords(latLonHPts,latLonRefEllips,aEllips,fEllips,is2D)
%%ELLIPS2OSCULATINGCOORDS Given points in [latitude;longitude;height]
%           ellipsoidal coordinates, convert it into the
%           [latitude;longitude;height] system of an osculating sphere
%           defined at the point latLonRefEllips on the surface of the
%           reference ellipsoid.
%
%INPUTS: latLonHPts The 3XN set of [latitude;longitude;height] points to
%               convert from ellipsoidal to osculating spherical
%               coordinates. Latitude and longitude are given in radians.
%               If a 2XN matrix is passed, it is assumed that the omitted
%               height coordinate is zero and the height is omitted from
%               latLonHOsc.
%   latLonRefEllips The [latitude;longitude] in radians in ellipsoidal
%               coordinates of where the osculating sphere is defined. This
%               is the point where the [latitude;longitude;height] in both
%               coordinate systems is the same (height=0).
%             a The semi-major axis length of the ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84SemiMajorAxis is used.
%             f The flattening factor of the ellipsoid. If this argument is
%               omitted or an empty matrix is passed, the value in
%               Constants.WGS84Flattening is used.
%          is2D Indicates whether one just cares about a 2D cut of the
%               ellipsoid. If so, then the radius used will be the radius
%               of curvature in the meridian. The default if omitted or an
%               empty matrix is passed is false.
%
%OUTPUTS: latLonHOsc A 3XN set of points in [latitude;longitude;height] in
%                    the osculating sphere coordinate system. If latLonHPts
%                    was 2XN, then the height component of this matrix is
%                    also omitted.
%       r, spherCent The range of the opsculating sphere and the center of
%                    the sphere given in the Cartesian coordinates
%                    associated with the reference ellipse.
%
%EXAMPLE:
%In this example, we show that the ellips2OsculatingCoords and
%osculatingCoords2Ellipse conversions are consistent. That is, converting
%to the osculating ellipse coordinate system and then converting back
%produces results that are within finite precision limits of the original
%point.
% latLonRef=deg2rad([20.776659;-156.010806]);
% latLonHPt=[deg2rad([20.267350;-155.837931]);10];
% latLonHOsc=ellips2OsculatingCoords(latLonHPt,latLonRef);
% RelErr=max(abs((osculatingCoords2Ellipse(latLonHOsc,latLonRef)-latLonHPt)./latLonHPt))
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<5||isempty(is2D))
        is2D=false;
    end
    
    if(nargin<4||isempty(fEllips))
        fEllips=Constants.WGS84Flattening;
    end

    if(nargin<3||isempty(aEllips))
        aEllips=Constants.WGS84SemiMajorAxis;
    end

    [r,spherCent]=osculatingSpher4LatLon(latLonRefEllips,aEllips,fEllips,is2D);
    ptsCartEllips=ellips2Cart(latLonHPts,aEllips,fEllips);
    ptsCartSpher=bsxfun(@minus,ptsCartEllips,spherCent);
    latLonHOsc=Cart2Ellipse(ptsCartSpher,[],r,0);
    
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
