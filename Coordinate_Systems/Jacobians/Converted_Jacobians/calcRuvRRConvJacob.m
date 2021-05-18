function J=calcRuvRRConvJacob(zRUVRR,useHalfRange,lTx,lRx,M)
%%CALCRUVRRCONVJACOB Calculate the Jacobian for a monostatic or bistatic
%            range direction cosines measurement in 3D with range rate,
%            ignoring atmospheric effects, with respect to Cartesian
%            position augmented with a range rate component. This type of
%            Jacobian is useful when performing tracking using Cartesian-
%            converted measurements where the clutter density is specified
%            in the measurement coordinate system, not the converted
%            measurement coordinate system.
%
%INPUTS: zRUVRR A 4X1 point in bistatic range, u and v, and range rate
%           coordinates in the format [bistatic range;u;v; range rate]. It
%           is also possible to pass a 3X1 point of just
%           [bistatic range;u;v], because the Jacobian does not depend on
%           the range-rate value.
% useHalfRange A boolean va0lue specifying whether the bistatic range value
%           should be divided by two. This normally comes up when operating
%           in monostatic mode, so that the range reported is a one-way
%           range. The default if this parameter is not provided, or an
%           empty matrix is passed, is false.
%       lTx The 3X1 [x;y;z] location vector of the transmitter in global
%           Cartesian coordinates. If this parameter is omitted or an
%           empty matrix is passed, then the transmitter is assumed to be
%           at the origin.
%       lRx The 3X1 [x;y;z] location vector of the receiver in Cartesian
%           coordinates. If this parameter is omitted or an empty matrix
%           is passed, then the receiver is assumed to be at the origin.
%         M A 3X3 rotation matrix to go from the alignment of the global
%           coordinate system to that at the receiver. The z-axis of the
%           local coordinate system of the receiver is the pointing
%           direction of the receiver. If omitted or an empty matrix is
%           passed, then it is assumed that the local coordinate system is
%           aligned with the global and M=eye(3) --the identity matrix is
%           used. 
%
%OUTPUTS: J The 4X4 Jacobian matrix. Each row is a components of
%           [range;u;v;range rate] in that order with derivatives taken
%           with respect to [x,y,z,range rate] by column.
%
%This function calls the calcRuvConvJacob function and augments the
%returned matrix with a 1 for the range-rate-range-rate derivative.
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

J=[calcRuvConvJacob(zRUVRR(1:3),useHalfRange,lTx(1:3),lRx(1:3),M),zeros(3,1);
    zeros(1,3),                                                   1];

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
