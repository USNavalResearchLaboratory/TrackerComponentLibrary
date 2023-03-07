function J=calcPolarJacob(x,systemType,useHalfRange,lTx,lRx)
%%CALCPOLARJACOB Calculate the Jacobian for a 2D polar measurement,
%                ignoring atmospheric effects with respect to Cartesian
%                position.
%
%INPUTS: x The 2XnumPoints positions of the target in Cartesian coordinates
%          in the order [x;y] where Jacobians are desired.
% systemType An optional parameter specifying the axis from which the
%          angles are measured. Possible values are
%          0 (The default if omitted or an empty matrix is passed) The
%            azimuth angle is counterclockwise from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided, or an
%          empty matrix is passed, is true.
%      lTx The 2X1 [x;y] location vector of the transmitter in global
%          Cartesian coordinates. If this parameter is omitted or an
%          empty matrix is passed, then the transmitter is assumed to be at
%          the origin.
%      lRx The 2X1 [x;y] location vector of the receiver in global
%          Cartesian coordinates. If this parameter is omitted or an
%          empty matrix is passed, then the receiver is assumed to be at
%          the origin.
%
%OUTPUTS: J The 2X2XnumPoints set of Jacobian matrices, one for each point
%           in x, with derivatives with respect to position components.
%           Each row is a component of range and azimuth in that order with
%           derivatives taken with respect to [x,y] across columns.
%
%This function just calls rangeGradient and polAngGradient. A rotation of
%the angular coordinate system by M does not affect the moments, so no
%rotation matrix is needed
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(lRx))
	lRx=zeros(2,1); 
end

if(nargin<4||isempty(lTx))
    lTx=zeros(2,1);
end

if(nargin<3||isempty(useHalfRange))
	useHalfRange=true; 
end

if(nargin<2||isempty(systemType))
	systemType=0; 
end

J=[rangeGradient(x(1:2,:),useHalfRange,lTx(1:2),lRx(1:2));polAngGradient(x(1:2,:),systemType,lRx(1:2))];

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
