function J=TDOAGradient(x,lRx1,lRx2,c)
%%TDOAGRADIENT Determine the gradient of a time-delay of arrival (TDOA)
%           measurement in 2D or 3D space with respect to position
%           (gradient components for velocity etc. are zero and are not
%           provided). Atmospheric and other propagation effects are not
%           taken into account.
%
%INPUTS: x The numPosDimX1 target position vector of the form [x;y] or
%          [x;y;z].
% lRx1,lRx2 The numPosDimX1 locations of the receivers. This assumes that
%          the TDOA measurement is taken using the time received at lRx1
%          as the reference time. Thus, the TDOA is the time at lRx1 minus
%          the time at lRx2. If either of these are omitted or an empty
%          matrix is passed, then the receiver is placed at the origin.
%          These cannot be both placed in the same spot.
%        c The propagation speed in the medium in question. If this
%          parameter is omitted or an empty matrix is passed, the default
%          value of Constants.speedOfLight is used.
%
%OUTPUTS: J A 1XnumPosDim gradient of the TDOA with derivatives taken with
%           respect to components [x,y,z] in 3D or [x,y] in 2D in that
%           order.
%
%A TDOA measurement is of the form (norm(x-lRx1)-norm(x-lRx2))*(1/c). The
%gradient is straightforward to find from there.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(x,1);

if(nargin<4||isempty(c))
   c=Constants.speedOfLight;
end

if(nargin<3||isempty(lRx2))
    lRx2=zeros(numDim,1); 
end

if(nargin<2||isempty(lRx1))
    lRx1=zeros(numDim,1); 
end

deltaRx1=x-lRx1;
deltaRx2=x-lRx2;

J=(1/c)*(deltaRx1.'/norm(deltaRx1)+deltaRx2.'/norm(deltaRx2));

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
