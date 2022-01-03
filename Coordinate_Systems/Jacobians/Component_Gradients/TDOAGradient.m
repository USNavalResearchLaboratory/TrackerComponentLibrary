function J=TDOAGradient(x,lRef,lRx,c)
%%TDOAGRADIENT Determine the gradient of a time-delay of arrival (TDOA)
%           measurement in 1D, 2D, or 3D space with respect to position
%           (gradient components for velocity etc. are zero and are not
%           provided). Atmospheric and other propagation effects are not
%           taken into account.
%
%INPUTS: x The numPosDimX1 target position vector of the form [x], [x;y] or
%          [x;y;z].
% lRef,lRx The numPosDimX1 locations of the receivers. This assumes that
%          the TDOA measurement is taken using the time received at lRx1
%          as the reference time. Thus, the TDOA is the time at lRx2 minus
%          the time at lRef. If either of these are omitted or an empty
%          matrix is passed, then the receiver is placed at the origin.
%          These cannot be both placed in the same spot.
%        c The propagation speed in the medium in question. If this
%          parameter is omitted or an empty matrix is passed, the default
%          value of Constants.speedOfLight is used.
%
%OUTPUTS: J A 1XnumPosDim gradient of the TDOA with partial derivatives
%           taken with respect to components [x,y,z] in 3D, [x,y] in 2D or
%           [x] in 1D, in that order.
%
%A TDOA measurement is of the form (norm(x-lRx1)-norm(x-lRx2))*(1/c). The
%gradient is straightforward to find from there.
%
%EXAMPLE:
%Here, we validate that the gradient is consistent with the value obtained
%by numerical differentiation.
% lRef=[-3;0];
% lRx=[3;0];
% c=1;
% xTrue=[4;4];
% h=@(x)getTDOA(x,lRef,lRx,c);
% grad=TDOAGradient(xTrue,lRef,lRx,c);
% gradNumDiff=numDiff(xTrue,h,1);
% relErr=(grad-gradNumDiff)./abs(grad)
%One will see that the relative error is better than 1e-10 in each
%component, indicating what that results are consistent.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(x,1);

if(nargin<4||isempty(c))
    c=Constants.speedOfLight;
end

if(nargin<3||isempty(lRx))
    lRx=zeros(numDim,1); 
end

if(nargin<2||isempty(lRef))
    lRef=zeros(numDim,1); 
end

deltaRx1=x-lRef;
deltaRx2=x-lRx;

J=(1/c)*(deltaRx2.'/norm(deltaRx2)-deltaRx1.'/norm(deltaRx1));

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
