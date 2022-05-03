function J=calcRuvRRJacob(x,useHalfRange,xTx,xRx,M)
%%CALCRUVRRJACOB Calculate the Jacobian of a 3D bistatic r-u-v measurement
%                with non-relativistic range rate, ignoring atmospheric
%                effects. The derivatives are taken with respect to the
%                Cartesian position and velocity components of the state.
%
%INPUTS: x The 6X1 position and velocity vector of the target in Cartesian
%          coordinates in the order [x;y;z;xDot;yDot;zDot].
% useHalfRange A boolean value specifying whether the bistatic range value
%          (and thus the range rate) should be divided by two. This
%          normally comes up when operating in monostatic mode, so that the
%          range reported is a one-way range. The default if this parameter
%          is not provided is false.
%      xTx The 6X1 position and velocity vector of the transmitter. The
%          default if omitted or an empty matrix is passed is a vector of
%          all zeros.
%      xRx The 6X1 position and velocity vector of the receiver. The
%          default if omitted or an empty matrix is passed is a vector of
%          all zeros.
%        M A 3X3 rotation matrices to go from the alignment of the global
%          coordinate system to that at the receiver. If omitted or an
%          empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(3) --the
%          identity matrix is used.
%
%OUTPUTS: J The 4X6 Jacobian matrix with derivatives with respect to
%           position components and velocity. Each row is a component of
%           [range;u;v;range rate] in that order with derivatives taken
%           with respect to [x,y,z,xDot,yDot,zDot] across columns.
%
%This function just calls the functions rangeGradient, uvGradient, and
%rangeRateGradient.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(M))
	M=eye(3,3); 
end

if(nargin<4||isempty(xRx))
	xRx=zeros(6,1); 
end

if(nargin<3||isempty(xTx))
    xTx=zeros(6,1);
end

if(nargin<2||isempty(useHalfRange))
	useHalfRange=false; 
end

J=zeros(4,6);
J(1,1:3)=rangeGradient(x(1:3),useHalfRange,xTx(1:3),xRx(1:3));
J(2:3,1:3)=uvGradient(x(1:3),xRx(1:3),M);
J(4,:)=rangeRateGradient(x(1:6),useHalfRange,xTx,xRx);

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

