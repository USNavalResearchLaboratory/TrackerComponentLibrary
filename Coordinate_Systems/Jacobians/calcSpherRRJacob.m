function J=calcSpherRRJacob(x,systemType,useHalfRange,lTx,lRx,M)
%%CALCSPHERRRJACOB Calculate the Jacobian of a 3D bistatic spherical
%                measurement with non-relativistic range rate, ignoring
%                atmospheric effects. The derivatives are taken with
%                respect to Cartesian position and velocity components of
%                the state.
%
%INPUTS: x The 6X1 position and velocity vector of the target in Cartesian
%          coordinates in the order [x;y;z;xDot;yDot;zDot].
% systemType An optional parameter specifying the axis from which the
%          angles are measured in radians. Possible values are
%          0 (The default if omitted) Azimuth is measured 
%            counterclockwise from the x-axis in the x-y plane. Elevation
%            is measured up from the x-y plane (towards the z-axis). This
%            is consistent with common spherical coordinate systems for
%            specifying longitude (azimuth) and geocentric latitude
%            (elevation).
%          1 Azimuth is measured counterclockwise from the z-axis in the
%            z-x plane. Elevation is measured up from the z-x plane
%            (towards the y-axis). This is consistent with some spherical
%            coordinate systems that use the z axis as the boresight
%            direction of the radar.
%          2 This is the same as 0 except instead of being given elevation,
%            one desires the angle away from the z-axis, which is
%            (pi/2-elevation).
%          3 This is the same as 0 except azimuth is measured clockwise
%            from the y-axis in the x-y plane instead of counterclockwise
%            from the x-axis. This coordinate system often arises when
%            given "bearings" in a local East-North-Up coordinate system,
%            where the bearing directions are measured East of North.
% useHalfRange An optional boolean value specifying whether the bistatic
%          (round-trip) range value has been divided by two. This normally
%          comes up when operating in monostatic mode (the most common
%          type of spherical coordinate system), so that the range
%          reported is a one-way range (or just half a bistatic range).
%          The default if this parameter is not provided is false if lTx
%          is provided and true if it is omitted (monostatic). 
%      lTx The 6X1 position and velocity vector of the transmitter.
%      lRx The 6X1 position and velocity vector of the receiver.
%        M A 3X3 rotation matrices to go from the alignment of the global
%          coordinate system to that at the receiver. If omitted, then it
%          is assumed that the local coordinate system is aligned with the
%          global and M=eye(3) --the identity matrix is used.
%
%OUTPUTS: J A 4X6 Jacobian matrix where the rows are
%          [bistatic range;azimuth;elevation; range rate] in that order and
%          the columns take the derivative of the row component with
%          respect to [x,y,z,xDot,yDot,zDot] in that order.
%
%This function just calls rangeGradient, spherAngGradient, and
%rangeRateGradient.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(M))
   M=eye(3,3); 
end

if((nargin<4||isempty(lTx))&&(nargin<3||isempty(useHalfRange)))
    useHalfRange=true;
elseif(nargin<3||isempty(useHalfRange))
    useHalfRange=false;
end

if(nargin<5||isempty(lRx))
    lRx=zeros(6,1);
end

if(nargin<4||isempty(lTx))
    lTx=zeros(6,1);
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

J=zeros(4,6);
J(1,1:3)=rangeGradient(x(1:3),useHalfRange,lTx(1:3),lRx(1:3));
J(2:3,1:3)=spherAngGradient(x(1:3),systemType,lRx(1:3),M);
J(4,:)=rangeRateGradient(x,useHalfRange,lTx,lRx);
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
