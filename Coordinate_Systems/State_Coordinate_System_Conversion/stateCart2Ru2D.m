function ru2DStates=stateCart2Ru2D(x)
%%STATECART2RU2D Convert a target state in 2D space consiting of position,
%               velocity and possibly acceleration terms into a monostatic
%               range and a single direction cosine coordines (assuming the
%               object is in front of the sensor). With acceleration, this
%               would be going from [x;y;xDot;yDot;xDDot;yDDot] to
%               [r;u;rDot;uDot;rDDot;uDot] where two Ds indicate a second
%               derivative with respect to time.
%
%INPUTS: x The 4XN or 6XN set of polar state vectors consisting of
%          position, velocity and possibly acceleration.
%
%OUTPUTS: ru2DStates The 4XN or 6XN set of target states given in polar
%                   coordinates. The range is a one-way range.
%
%The derivation is given in [1]. The function aCVRu2D is an implementation
%of a related linear dynamic model.
%
%EXAMPLE:
%Here we note that the results are consistent with the inverse function:
%stateRu2D2Cart.
% x=[100;60;3;12;108;-116];
% xRet=stateRu2D2Cart(stateCart2Ru2D(x));
% max(abs(x(:)-xRet(:)))
%One will see that the error is less than 1e-13, indicating good agreement.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic Linear Cartesian Dynamic Models in Local
%    Coordinates," Naval Research Laboratory7: Washignton, DC, No.
%    NRL/MR/5344--19-9882, 24 Aug. 2019.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

ruPos=Cart2Ru2D(x(1:2,:),true);

numDim=size(x,1);
N=size(x,2);

ru2DStates=zeros(numDim,N);

ru2DStates(1:2,:)=ruPos;

r=ruPos(1,:);
u=ruPos(2,:);

v2=max(0,1-u.^2);
v=sqrt(v2);
uVecr=[u;v];
uVecu=[v;-u];

rVel=x(3:4,:);
rDot=sum(rVel.*uVecr,1);
uDot=v.*sum(rVel.*uVecu,1)./r;

ru2DStates(3:4,:)=[rDot;uDot];

if(numDim>4)
    rAccel=x(5:6,:);
    ar=sum(rAccel.*uVecr,1);
    au=sum(rAccel.*uVecu,1);
    
    rDDot=ar+r.*uDot.^2./v2;
    uDDot=(1./r).*(au.*v-2*rDot.*uDot)-u.*uDot.^2./v2;
    
    ru2DStates(5:6,:)=[rDDot;uDDot];
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
