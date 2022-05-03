function z=Cart2TDOASphere(points,systemType,lRef,lRx,M,c)
%%CART2TDOASPHERE Given Cartesion positions, convert them to a
%       time-difference of arrival (TDOA) component and spherical angles
%       The angles are assumed taken at the non-reference receiver.
%
%INPUTS: points A 3XN set of Cartesian points.
% systemType An optional parameter specifying the axis from which the
%           angles are measured in radians. Possible values are
%           0 (The default if omitted) Azimuth is measured 
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z-axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one is given the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
%      lRef The 3X1 location of the reference receiver in global
%           coordinates.
%       lRx The 3X1 location of the receiver in global coordinates at which
%           the direction of arrival in spherical coordinates is locally
%           defined and whose time of arrival (TOA) has the TOA at lRef
%           subtracted to get the TDOA.
%         M A 3X3 rotation matrix to go from the alignment of the global
%           coordinate system to that at the receiver. If omitted, then it
%           is assumed that the local coordinate system is aligned with the
%           global and M=eye(3) --the identity matrix is used.
%         c The propagation speed in the medium in question. If this
%           parameter is omitted or an empty matrix is passed, the default
%           value of Constants.speedOfLight is used.
%
%OUTPUTS: z A 3XN set of points of the form [TDOA;azimuth;elevation] with
%           the angles given in radians.
%
%This function just calls getTDOA and getSpherAngle.
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(c))
    c=Constants.speedOfLight;
end

if(nargin<5||isempty(M))
    M=eye(3,3);
end

if(isempty(systemType))
    systemType=0;
end

z=[getTDOA(points,lRef,lRx,c);
   getSpherAngle(points,systemType,lRx,M)];

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
