function v=headingSpeed2Vel(plhPoint,speed,geoEastOfNorth,angUpFromLevel)
%%HEADINGSPEED2VEL Given a target location in latitude and logitude a
%    target speed, the heading East of North and any angle it is making up
%    from level flight, obtain a 3D velocity vector. This function can be
%    useful for parameterizing 
%
%INPUTS: plhPoint The 2XN or 3XN location matrix of the points in geodetic
%              latitude  and longitude in radians at which the headings are
%              taken. The point can be [latitude;longitude] or
%              [latitude;longitude;height]. The height component is ignored
%              if included because it does not change the result.
%        speed An NX1 or 1XN array of the speed of the target associated
%              with each point. If a single scalar is passed, then the
%              speed is assumed to be the same for all targets.
% geoEastOfNorth An NX1 or 1XN array of N geographic headings in radians
%              clockwise from North that should be turned into ECEF unit
%              vectors. A geographic heading is a direction in the local
%              tangent plane of an East-North-Up coordinate system as
%              defined on a specific reference ellipsoid. If all headings
%              are the same (but the angles up from level vary), then this
%              can be a scalar value. If this parameter is omitted or an
%              empty matrix is passed, then the default of 0 is used.
% angUpFromLevel An NX1 or 1XN array of N angles of the trajectory above
%              the local tangent plane to the reference ellipsoid. If all
%              elevations are the same (but geographic headings might vary),
%              then this can be a scalar value. If this parameter is
%              omitted or an empty matrix is passed, then the default of 0
%              is used.
%
%OUTPUTS: v A 3XN matrix of velocity vectors. The units are the same as the
%           units of the speed input.
%
%This function just calls geogHeading2uVec and multiplies the result by
%speed.
%
%May 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(angUpFromLevel))
    angUpFromLevel=0; 
end

if(nargin<3||isempty(geoEastOfNorth))
    geoEastOfNorth=0; 
end

v=bsxfun(@times,speed(:).',geogHeading2uVec(plhPoint,geoEastOfNorth,angUpFromLevel));

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
