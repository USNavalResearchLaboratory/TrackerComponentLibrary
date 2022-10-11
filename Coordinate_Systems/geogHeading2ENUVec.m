function u=geogHeading2ENUVec(heading,angUpFromLevel)
%%GEOGHEADING2ENUVEC Given a geographic heading in radians East of North
%       and an elevation above the local tangent plane, get a corresponding
%       unit vector in the specified direction in 3D ENU coordinates.
%
%INPUTS: heading A length numPts vector of headings in radians East of
%                North.
% angUpFromLevel A length numPts vector of elevations above the local
%                tangent plane in radians.
%
%OUTPUTS: u A 3XnumPts set of unit vectors in local ENU coordinates
%           corresponding to the directions given in heading and
%           angUpFromLevel.
%
%The transformation is just a specific definition of the spherical
%coordinate system. There is no need for information on the shape of the
%reference ellipsoid.
%
%EXAMPLE:
%Here, we show that the values produced by this function are consistent
%with those returned by the ENUVec2GeogHeading function. The difference
%between the original heading and elevation and converted back values are
%on the order of finite precision errors.
% headingInit=[deg2rad(20);deg2rad(21)];
% angUpFromLevelInit=[deg2rad(5);deg2rad(51)];
% [heading,angUpFromLevel]=ENUVec2GeogHeading(geogHeading2ENUVec(headingInit,angUpFromLevelInit));
% abs(heading(:)-headingInit(:))
% abs(angUpFromLevel(:)-angUpFromLevelInit(:))
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

angUpFromLevel=angUpFromLevel(:).';
heading=heading(:).';

cosEl=cos(angUpFromLevel);
sinEl=sin(angUpFromLevel);

u=[sin(heading).*cosEl;
   cos(heading).*cosEl;
   sinEl];

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
