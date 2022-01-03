function TDOA=getTDOA(x,lRef,lRx,c)
%%FUNCTIONGETTDOA Determine the time-difference of arrival (TDOA) in 2D or
%           3D space with respect to position. Atmospheric and other
%           propagation effects are not taken into account.
%
%INPUTS: x The numPosDimXN set of target position vectors of the form
%          [x;y] or [x;y;z].
% lRef,lRx The numPosDimXN locations of the receivers. This assumes that
%          the TDOA measurement is taken using the time received at lRef
%          as the reference time. Thus, the TDOA is the time at lRx minus
%          the time at lRef. If either of these are omitted or an empty
%          matrix is passed, then the receivers are placed at the origin.
%          These cannot be both placed in the same spot. If a receiver is
%          the same for all N positions in x, then a single numPosDimX1
%          vector can be passed for its position.
%        c The propagation speed in the medium in question. If this
%          parameter is omitted or an empty matrix is passed, the default
%          value of Constants.speedOfLight is used.
%
%OUTPUTS: TDOA A 1XN vector contianing all of the TDOA values.
%
%TDOA measurements commonly arise when performing passive localization.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(x,2);
numDim=size(x,1);

if(nargin<4||isempty(c))
   c=Constants.speedOfLight;
end

if(nargin<3||isempty(lRx))
    lRx=zeros(numDim,N);
elseif(size(lRx,2)==1)
    lRx=repmat(lRx,[1,N]);
end

if(nargin<2||isempty(lRef))
    lRef=zeros(numDim,N);
elseif(size(lRef,2)==1)
    lRef=repmat(lRef,[1,N]);
end

TDOA=(1/c)*(sqrt(sum((x-lRx).^2,1))-sqrt(sum((x-lRef).^2,1)));
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
