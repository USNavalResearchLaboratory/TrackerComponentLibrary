function J=calcSpherConvJacob(zSpher,systemType,useHalfRange,lTx,lRx,M)
%%CALCSPHERCONVJACOB Calculate the Jacobian for a monostatic or bistatic
%            range and a spherical direction measurement in 3D, ignoring
%            atmospheric effects, with respect to Cartesian position. This
%            type of Jacobian is useful when performing tracking using
%            Cartesian-converted measurements where the clutter density is
%            specified in the measurement coordinate system, not the
%            converted measurement coordinate system.
%
%INPUTS: zSpher A 3XN set of points in range and azimuth and angle in the
%           format [range;azimuth;angle], where the angles are given in
%           radians.
% systemType An optional parameter specifying the axes from which the
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
%             coordinate systems that use the z axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one desires the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
% useHalfRange An optional boolean value specifying whether the bistatic
%           (round-trip) range value has been divided by two. This normally
%           comes up when operating in monostatic mode (the most common
%           type of spherical coordinate system), so that the range
%           reported is a one-way range (or just half a bistatic range).
%           The default if this parameter is not provided is false if lTx
%           is provided and true if it is are omitted (monostatic).
%       lTx The 3X1 [x;y;z] location vector of the transmitter in global
%           Cartesian coordinates. If this parameter is omitted or an
%           empty matrix is passed, then the transmitter is assumed to be
%           at the origin.
%       lRx The 3X1 [x;y;z] location vector of the receiver in Cartesian
%           coordinates. If this parameter is omitted or an empty matrix
%           is passed, then the receiver is assumed to be at the origin.
%         M A 3X3 rotation matrices to go from the alignment of the global
%           coordinate system to that at the receiver. If omitted, then it
%           is assumed that the local coordinate system is aligned with the
%           global and M=eye(3) --the identity matrix is used.
%
%OUTPUTS: J The 3X3XN set of Jacobian matrices, one for each point given.
%           Each row is a components of [range;azimuth;elevation] in that
%           order with derivatives taken with respect to [x,y,z] by column.
%
%This function converts the measurement into Cartesian coordinates and then
%calls rangeGradient and spherAngGradient. Note that singularities exist at
%the poles.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(M))
	M=eye(3,3); 
end

if(nargin<5)
   lRx=[]; 
end

if(nargin<4)
    lTx=[];
end

if(isempty(lTx)&&(nargin<3||isempty(useHalfRange)))
    useHalfRange=true;
elseif(~isempty(lTx)&&(nargin<3||isempty(useHalfRange)))
    useHalfRange=false;
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

N=size(zSpher,2);

x=spher2Cart(zSpher,systemType,useHalfRange,lTx,lRx,M);

J=zeros(3,3,N);
J(1,:,:)=rangeGradient(x,useHalfRange,lTx,lRx);
J(2:3,:,:)=spherAngGradient(x,systemType,lRx,M);

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
