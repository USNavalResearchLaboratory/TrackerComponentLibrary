function H=calcSpherHessian(x,systemType,useHalfRange,lTx,lRx,M)
%%CALCSPHERHESSIAN Calculate the Hessian matrix (a matrix of second partial
%          derivatives) for a monostatic or bistatic spherical measurement,
%          ignoring atmospheric effects, with respect to 3D Cartesian
%          position.
%
%INPUTS: x The 3X1 position of the target in Cartesian coordinates in the
%          order [x;y;z].
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
%          2 This is the same as 0 except instead of being given
%            elevation, one desires the angle away from the z-axis, which
%            is (pi/2-elevation).
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
%      lTx The 3X1 transmitter position in the global coordinate system
%          with [x;y;z] components. If omitted or an empty matrix is
%          passed, then a vector of zeros is used.
%      lRx The 3X1 receiver position in the global coordinate system
%          with [x;y;z] components. If omitted or an empty matrix is
%          passed, then a vector of zeros is used.
%        M A rotation matrix from the global Coordinate system to the
%          orientation of the coordinate system at the receiver. This is
%          only necessary if UV direction components are desired. If
%          omitted, it is assumed to be the identity matrix.
%
%OUTPUTS: H The 3X3X3 Hessian matrix, where H(:,:,1) is the Hessian for the
%           range component, H(:,:,2) is the Hessian for the azimuth
%           component, and H(:,:,3) is the Hessian for the elevation
%           component. The ordering of the derivatives in each matrix is:
%                  [d^2/(dxdx), d^2/(dxdy), d^2/(dxdz);
%                   d^2/(dydx), d^2/(dydy), d^2/(dydz);
%                   d^2/(dzdx), d^2/(dzdy), d^2/(dzdz)];
%          note that each matrix is symmetric (i.e.
%                   d^2/(dydx)=d^2/(dxdy) ).
%
%This function just calls rangeHessian and spherAngHessian.
%
%June 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
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
    lRx=zeros(3,1);
end

if(nargin<4||isempty(lTx))
    lTx=zeros(3,1);
end

if(nargin<2||isempty(systemType))
    systemType=0;
end

H=zeros(3,3,3);
H(:,:,1)=rangeHessian(x,useHalfRange,lTx,lRx);
H(:,:,2:3)=spherAngHessian(x,systemType,lRx,M);
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
