function J=calcRuvJacob(x,useHalfRange,lTx,lRx,M)
%%CALCRUVJACOB Calculate the Jacobian for a monostatic or bistatic r-u-v
%           measurement with respect to 3D Cartesian position. Atmospheric
%           effects are ignored.
%
%INPUTS: x The 3X1 position of the target in Cartesian coordinates in the
%          order [x;y;z].
% useHalfRange A boolean value specifying whether the bistatic range value
%          has been divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided is false.
%      lTx The 3X1 transmitter position in the global coordinate system
%          with [x;y;z] components. If omitted or an empty matrix is
%          passed, then a vector of zeros is used.
%      lRx The 3X1 receiver position in the global coordinate system
%          with [x;y;z] components. If omitted or an empty matrix is
%          passed, then a vector of zeros is used.
%        M A 3X3 rotation matrix from the global Coordinate system to the
%          orientation of the coordinate system at the receiver. This is
%          only necessary if u-v direction components are desired. If
%          omitted, it is assumed to be the identity matrix.
%
%OUTPUTS: J The 3X3 Jacobian matrix with derivatives with respect to
%           position components. Each row is a component of bistatic range,
%           u and v in that order with derivatives taken with respect to
%           [x,y,z] across columns.
%
%This function just calls rangeGradient and uvGradient.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(M))
	M=eye(3,3); 
end

if(nargin<4||isempty(lRx))
	lRx=zeros(3,1); 
end

if(nargin<3||isempty(lTx))
    lTx=zeros(3,1);
end

if(nargin<2||isempty(useHalfRange))
	useHalfRange=false; 
end

J=zeros(3,3);
J(1,:)=rangeGradient(x(1:3),useHalfRange,lTx(1:3),lRx(1:3));
J(2:3,:)=uvGradient(x(1:3),lRx(1:3),M);

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

