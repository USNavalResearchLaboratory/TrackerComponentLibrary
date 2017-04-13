function J=rangeGradient(x,useHalfRange,lTx,lRx)
%%RANGEGRADIENT Determine the gradient of a 2D or 3D bistatic range
%           measurement with respect to position (gradient components for
%           velocity etc. are zero and are not provided). Atmospheric and
%           other propagation effects are not taken into account.
%
%INPUTS: x A numPosDimX1 target position vector of the form [x;y] or
%          [x;y;z].
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up
%          when operating in monostatic mode, so that the range reported is
%          a one-way range. The default if this parameter is not provided
%          is false.
%      lTx The 3X1 (in 3D) or 2X1 (in 2D) position vector of the
%          transmitter. If this parameter is omitted or an empty matrix is
%          passed, then the transmitter is assumed to be at the origin.
%      lRx The 3X1 (in 3D) or 2X1 (in 2D) position vector of the receiver.
%          If this parameter is omitted or an empty matrix is passed, then
%          the receiver is assumed to be at the origin.
%
%OUTPUTS: J A 1XnumPosDim gradient of the bistatic range with derivatives
%           taken with respect to components [x,y,z] in 3D or [x,y] in 2D
%           in that order.
%
%Gradients with respect to bistatic range are discussed in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(x,1);

if(nargin<2||isempty(useHalfRange))
   useHalfRange=false; 
end

if(nargin<4||isempty(lRx))
    lRx=zeros(numDim,1); 
end

if(nargin<3||isempty(lTx))
    lTx=zeros(numDim,1); 
end

deltaRx=x-lRx;
deltaTx=x-lTx;

J=deltaRx.'/norm(deltaRx)+deltaTx.'/norm(deltaTx);

if(useHalfRange)
    J=J/2; 
end
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
