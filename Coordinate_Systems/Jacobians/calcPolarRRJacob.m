function J=calcPolarRRJacob(x,systemType,useHalfRange,lTx,lRx)
%%CALCPOLARRRJACOB Calculate the Jacobian of a 2D polar measurement,
%                  with non-relativistic range rate, ignoring atmospheric
%                  effects. The derivatives are taken with respect to
%                  Cartesian position and velocity components of the state.
%
%INPUTS: x The 4X1 position and velocity vector of the target in Cartesian
%          coordinates in the order [x;y;xDot;yDot].
% systemType An optional parameter specifying the axis from which the
%          angles are measured. Possible values are
%          0 (The default if omitted or an empty matrix is passed) The
%            azimuth angle is counterclockwise from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
% useHalfRange A boolean value specifying whether the bistatic range value
%          (and thus the range rate) should be divided by two. This
%          normally comes up when operating in monostatic mode, so that the
%          range reported is a one-way range. The default if this parameter
%          is not provided is false.
%      lTx The 4X1 position and velocity vector of the transmitter.
%      lRx The 4X1 position and velocity vector of the receiver.
%
%OUTPUTS: J The 2X4 Jacobian matrix with derivatives with respect to
%           position components and velocity. Each row is a component of
%           [range;azimuth;range rate] in that order with derivatives taken
%           with respect to [x,y,xDot,yDot] across columns.
%
%This function just calls the functions rangeGradient, polAngGradient and
%rangeRateGradient. The rotation of the angular reference point does not
%matter (i.e. the definition of the local x direction).
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(lRx))
	lRx=zeros(4,1); 
end

if(nargin<4||isempty(lTx))
    lTx=zeros(4,1);
end

if(nargin<3||isempty(useHalfRange))
	useHalfRange=true; 
end

if(nargin<2||isempty(systemType))
	systemType=0; 
end

J=zeros(2,4);
J(1,1:2)=rangeGradient(x(1:2),useHalfRange,lTx(1:2),lRx(1:2));
J(2,1:2)=polAngGradient(x(1:2),systemType,lRx(1:2));
J(3,:)=rangeRateGradient(x,useHalfRange,lTx,lRx);
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
