function polStates=stateCart2Pol(x,systemType)
%%STATECART2POL Convert a target state in 2D space consiting of position,
%               velocity and possibly acceleration terms into monostatic
%               polar coordinates. With acceleration, this would be going
%               from [x;y;xDot;yDot;xDDot;yDDot] to
%               [r;theta;rDot;thetaDot;rDDot;thetaDDot] where two Ds
%               indicates a second derivative with respect to time.
%
%INPUTS: x The 4XN or 6XN set of Cartesian state vectors consisting of
%          position, velocity and possibly acceleration. Note that a
%          singularity exists at the origin, so non-finite terms can arise
%          in the output for entry i where x(1:2,i)=[0;0].
% systemType An optional parameter specifying the axis from which the
%          angles are measured. Possible values are
%          0 (The default if omitted) The azimuth angle is counterclockwise
%            from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
%
%OUTPUTS: polStates The 4XN or 6XN set of target states given in polar
%                   coordinates. The range is a one-way range and the
%                   angles are given in radians.
%
%The derivation is given in [1]. The function aCVPolar is an implementation
%of a related linear dynamic model.
%
%EXAMPLE:
%Here we note that the results are consistent with the inverse function:
%statePol2Cart.
% systemType=0;
% x=[[100;-60;-3;12;108;-116],[1;1;1;1;1;1]];
% xRet=statePol2Cart(stateCart2Pol(x,systemType),systemType);
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

if(nargin<2||isempty(systemType))
    systemType=0; 
end

posPol=Cart2Pol(x(1:2,:),systemType,true);

numDim=size(x,1);
N=size(x,2);

polStates=zeros(numDim,N);

polStates(1:2,:)=posPol;

r=posPol(1,:);
ur=bsxfun(@rdivide,x(1:2,:),r);

if(systemType==0)
    cosTheta=ur(1,:);
    sinTheta=ur(2,:);
    
    uTheta=[-sinTheta;cosTheta];
else%systemType=1
    cosTheta=ur(2,:);
    sinTheta=ur(1,:);
    
    uTheta=[cosTheta;-sinTheta];
end

rVel=x(3:4,:);

rDot=sum(ur.*rVel,1);
thetaDot=sum(uTheta.*rVel,1)./r;

polStates(3,:)=rDot;
polStates(4,:)=thetaDot;

if(numDim>4)
    %If acceleration is provided.
    rAccel=x(5:6,:);
    ar=sum(ur.*rAccel,1);
    aTheta=sum(uTheta.*rAccel,1);
    
    rDDot=ar+r.*thetaDot.^2;
    thetaDDot=(1./r).*(aTheta-2*rDot.*thetaDot);
    
    polStates(5,:)=rDDot;
    polStates(6,:)=thetaDDot;
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
