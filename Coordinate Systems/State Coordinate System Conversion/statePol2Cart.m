function CartStates=statePol2Cart(x,systemType)
%%STATEPOL2CART Convert a target state in 2D space in monostatic polar
%               coordinates with first and possibly second time derivatives
%               of the position components into Cartesian coordinates. The
%               state has the format
%               [r;theta;rDot;thetaDot;rDDot;thetaDDot], where two Ds
%               indicate a second derivative with respect to time. In
%               Cartesian coordinates, the converted state has the form
%               [x;y;xDot;yDot;xDDot;yDDot].
%
%INPUTS: x The 4XN or 6XN set of polar state vectors consisting of
%          position, velocity and possibly acceleration. The angles are
%          given in radians. The range is a one-way range.
% systemType An optional parameter specifying the axis from which the
%          angles are measured. Possible values are
%          0 (The default if omitted) The azimuth angle is counterclockwise
%            from the x axis.
%          1 The azimuth angle is measured clockwise from the y axis.
%
%OUTPUTS: CartStates The 4XN or 6XN set of target states given in Cartesian
%                    coordinates.
%
%The derivation is given in [1]. The function aCVPolar is an implementation
%of a related linear dynamic model.
%
%EXAMPLE:
%Here we note that the results are consistent with the inverse function:
%stateCart2Pol.
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

numDim=size(x,1);
N=size(x,2);

CartStates=zeros(numDim,N);

r=x(1,:);
theta=x(2,:);
rDot=x(3,:);
thetaDot=x(4,:);

cosTheta=cos(theta);
sinTheta=sin(theta);

switch(systemType)
    case 0
        ur=[cosTheta;sinTheta];
        uTheta=[-sinTheta;cosTheta];
    case 1
        ur=[sinTheta;cosTheta];
        uTheta=[cosTheta;-sinTheta];
    otherwise
        error('Unknown system Type specified.')
end

%Position components.
CartStates(1:2,:)=bsxfun(@times,r,ur);
%Velocity components.
CartStates(3:4,:)=bsxfun(@times,rDot,ur)+bsxfun(@times,r.*thetaDot,uTheta);

if(numDim>4)
    %If acceleration is provided.
    rDDot=x(5,:);
    thetaDDot=x(6,:);
    
    CartStates(5:6,:)=bsxfun(@times,(rDDot-r.*thetaDot.^2),ur)+bsxfun(@times,(r.*thetaDDot+2*rDot.*thetaDot),uTheta);
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
