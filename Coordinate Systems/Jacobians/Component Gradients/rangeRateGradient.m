function J=rangeRateGradient(x,useHalfRange,lTx,lRx)
%%RANGERATEGRADIENT Determine the gradient of a 2D or 3D bistatic range
%          rate measurement with respect to the position and velocity
%          components of a target state. Higher order gradient terms are
%          not provided and are zero. Relativity and atmospheric effects
%          are not taken into account.
%
%INPUTS: x The 4X1 (for 2D position) or 6X1 (for 3D position) target state
%          with components of [x;y;z;xDot;yDot;zDot] in 3D or
%          [x;y;xDot;yDot] in 2D.
% useHalfRange A boolean value specifying whether the bistatic range value
%          (and thus the range rate) should be divided by two. This
%          normally comes up when operating in monostatic mode, so that the
%          range reported is a one-way range. The default if this parameter
%          is not provided is false.
%      lTx The 6X1 (in 3D) or 4X1 (in 2D) position and velocity vector of
%          the transmitter.
%      lRx The 6X1 (in 3D) or 4X1 (in 2D) position and velocity vector of
%          the receiver.
%
%OUTPUTS: J A 1X6 (in 3D) or 1X4 (in 2D) gradient of the bistatic range
%           rate with derivatives taken with respect to components
%           [x,y,z,xDot,yDot,zDot] in 3D or [x,y,xDot,yDot] in 2D in that
%           order.
%
%The non-relativistic range-rate approximation is derived in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(x,1);
numDim=xDim/2;

if(nargin<2||isempty(useHalfRange))
   useHalfRange=true; 
end

if(nargin<3||isempty(lTx))
    lTx=zeros(xDim,1); 
end

if(nargin<4||isempty(lRx))
    lRx=zeros(xDim,1); 
end

switch(numDim)
    case 2
        J=zeros(1,4);
    
        dtr=x(1:2)-lRx(1:2);
        dtl=x(1:2)-lTx(1:2);
        
        drB=dtr/norm(dtr)+dtl/norm(dtl);
        
        dvr=x(3:4)-lRx(3:4);
        dvl=x(3:4)-lTx(3:4);

        %The position derivaties
        J(1:2)=A2(dtr)*dvr/norm(dtr)^3+A2(dtl)*dvl/norm(dtl)^3;

        %The velocity derivatives.
        J(3:4)=drB;
    case 3
        J=zeros(1,6);
    
        dtr=x(1:3)-lRx(1:3);
        dtl=x(1:3)-lTx(1:3);
        
        drB=dtr/norm(dtr)+dtl/norm(dtl);
        
        dvr=x(4:6)-lRx(4:6);
        dvl=x(4:6)-lTx(4:6);

        %The position derivaties
        J(1:3)=A(dtr)*dvr/norm(dtr)^3+A(dtl)*dvl/norm(dtl)^3;

        %The velocity derivatives.
        J(4:6)=drB;
    otherwise
        error('Invalid state dimensionality.')
end

if(useHalfRange)
    J=J/2; 
end
end

function val=A(vec)
%%%A This is a simple helper function as defined in Appendix I of [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.

    x=vec(1);
    y=vec(2);
    z=vec(3);
    val=[y^2+z^2,-x*y,   -x*z;
        -x*y,   x^2+z^2,-y*z;
        -x*z,   -y*z,   x^2+y^2];
end

function val=A2(vec)
%%%A2 This is a 2D version of the helper function A as defined in Appendix
%     I of [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.

    x=vec(1);
    y=vec(2);
    val=[y^2, -x*y;
        -x*y,   x^2];
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
