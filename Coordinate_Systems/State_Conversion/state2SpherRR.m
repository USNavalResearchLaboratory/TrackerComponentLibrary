function z=state2SpherRR(xTar,systemType,useHalfRange,xTx,xRx,M)
%%STATE2SPHERRR Convert state vectors consisting of at least 3D position
%               and velocity in 3D space into local spherical coordinates
%               with non-relativistic range rate. An option allows for the
%               angles to be specified from different axes. Optionally, a
%               bistatic range can be used when considering bistatic
%               measurements in a local spherical coordinate system.
%
%INPUTS: xTar A xDimXN matrix of N target states consisting of 3D position
%             and velocity components in the order
%             xTar=[xPosition;xVelocity] and possible other components,
%             which will be ignored.
%  systemType An optional parameter specifying the axes from which
%             the angles are measured. Possible vaues are
%             0 (The default if omitted) Azimuth is measured
%               counterclockwise from the x-axis in the x-y plane.
%               Elevation is measured up from the x-y plane (towards the
%               z-axis). This is consistent with common spherical
%               coordinate systems for specifying longitude (azimuth)
%               and geocentric latitude (elevation).
%             1 Azimuth is measured counterclockwise from the z-axis in
%               the z-x plane. Elevation is measured up from the z-x plane
%               (towards the y-axis). This is consistent with some
%               spherical coordinate systems that use the z axis as
%               the boresight direction of the radar.
%             2 This is the same as 0 except instead of being given
%               elevation, one desires the angle away from the z-axis,
%               which is (pi/2-elevation).
%             3 This is the same as 0 except azimuth is measured clockwise
%               from the y-axis in the x-y plane instead of
%               counterclockwise from the x-axis. This coordinate system
%               often arises when given "bearings" in a local East-North-Up
%               coordinate system,  where the bearing directions are
%               measured East of North.
% useHalfRange An optional boolean value specifying whether the bistatic
%             (round-trip) range value has been divided by two. This
%             normally comes up when operating in monostatic mode (the most
%             common type of spherical coordinate system), so that the
%             range reported is a one-way range (or just half a bistatic
%             range). The default if this parameter is not provided is
%             false if xTx is provided and true if it is omitted
%             (monostatic).
%         xTx An xTxDimXN matrix of the states of the transmitters
%             consisting of stacked 3D position and velocity components.
%             Other components will be ignored. If this parameter is
%             omitted, the transmitters are assumed to be stationary at the
%             origin. If only a single vector is passed, then the
%             transmitter state is assumed the same for all of the target
%             states being converted.
%         xRx An xRxDimXN matrix of the states of the receivers
%             consisting of stacked 3D position and velocity components.
%             Other components will be ignored. If this parameter is
%             omitted, the receivers are assumed to be stationary at the
%             origin. If only a single vector is passed, then the
%             receiver state is assumed the same for all of the target
%             states being converted.
%           M A 3X3XN hypermatrix of the rotation matrices to go from the
%             alignment of the global coordinate system to that at the
%             receiver. If omitted, then it is assumed that the local
%             coordinate system is aligned with the global and M=eye(3).
%             If only a single 3X3 matrix is passed, then it is assumed to
%             be the same for all of the N conversions.
%
%OUTPUTS: z A 4XN matrix of the target states in xTar converted into
%           bistatic spherical and bistatic range rate coordinates. If
%           useHalfRange=true, then the r component is half the bistatic
%           range and the range rate is correspondingly halved.
%
%Details of the conversions are given in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(xTar,2);

if(nargin<6||isempty(M))
    M=repmat(eye(3),[1,1,N]);
elseif(size(M,3)==1&&N>1)
    M=repmat(M,[1,1,N]);
end

if(nargin<4)
    xTx=[];
end

if(nargin<5)
    xRx=[];
end

if(isempty(xTx)&&(nargin<3||isempty(useHalfRange)))
    useHalfRange=true;
elseif(nargin<3||isempty(useHalfRange))
    useHalfRange=false;
end

if(nargin<5||isempty(xRx))
    xRx=zeros(6,N);
elseif(size(xRx,2)==1&&N>1)
    xRx=repmat(xRx,[1,N]);
end

if(nargin<4||isempty(xTx))
    xTx=zeros(6,N);
elseif(size(xTx,2)==1&&N>1)
    xTx=repmat(xTx,[1,N]);
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

%Allocate space for the converted states.
z=zeros(4,N);

%Convert the positions.
z(1:3,:)=Cart2Sphere(xTar,systemType,useHalfRange,xTx,xRx,M);

%Compute the bistatic range rates.
z(4,:)=getRangeRate(xTar(1:6,:),useHalfRange,xTx,xRx);
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
