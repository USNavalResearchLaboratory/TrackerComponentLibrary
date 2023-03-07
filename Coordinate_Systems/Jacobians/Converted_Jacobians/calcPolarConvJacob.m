function J=calcPolarConvJacob(zPolar,systemType,useHalfRange,lTx,lRx,M)
%%CALCPOLARCONVJACOB Calculate the Jacobian for a monostatic or bistatic
%            range and polar angle measurement in 2D with respect to
%            Cartesian position. Atmospheric effects are ignored. This
%            type of Jacobian is useful when performing tracking using
%            Cartesian-converted measurements where the clutter density is
%            specified in the measurement coordinate system, not the
%            converted measurement coordinate system.
%
%INPUTS: zPolar A 2X1 point in polar coordinates in the format
%          [range;azimuth], where the angle is given in radians and the
%          range can be bistatic.
% systemType An optional parameter specifying the axis from which the
%          azimuth angle is measured. It is assumed that the azimuth
%          angle is given in radians. Possible values are
%          0 (The default if omitted) The azimuth angle is
%             counterclockwise from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up
%          when operating in monostatic mode, so that the range reported is
%          a one-way range. The default if this parameter is not provided
%          is false.
%      lTx The 2X1 [x;y] location vector of the transmitter in Cartesian
%          coordinates. If this parameter is omitted or an empty matrix is
%          passed, then the transmitter is assumed to be at the origin.
%      lRx The 2X1 [x;y] location vector of the receiver in Cartesian
%          coordinates. If this parameter is omitted or an empty matrix is
%          passed, then the receiver is assumed to be at the origin.
%        M A 2X2 rotation matrices to go from the alignment of the global
%          coordinate system to that at the receiver. If omitted or an
%          empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(2,2)
%          --the identity matrix is used. 
%
%OUTPUTS: J The 2X2 Jacobian matrix. Each row is a components of range, and
%           azimuth (in that order by row) with derivatives taken with
%           respect to [x,y] by column.
%
%This function converts the measurement into Cartesian coordinates and then
%calls rangeGradient and polAngGradient.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(M))
   M=eye(2,2); 
end

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

x=pol2Cart(zPolar,systemType,useHalfRange,lTx,lRx,M);

J=zeros(2,2);
J(1,:)=rangeGradient(x,useHalfRange,lTx,lRx);
J(2,:)=polAngGradient(x,systemType,lRx);

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
