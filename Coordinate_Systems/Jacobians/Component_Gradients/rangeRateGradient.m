function JTotal=rangeRateGradient(xG,useHalfRange,lTx,lRx)
%%RANGERATEGRADIENT Determine the gradient of a 1D, 2D or 3D bistatic range
%          rate measurement with respect to the position and velocity
%          components of a target state. Higher order gradient terms are
%          not provided and are zero. Relativity and atmospheric effects
%          are not taken into account. If the target and transmitter are
%          collocated, the target is assumed to be the transmitter.
%
%INPUTS: xG The (2*numDims)XN target state vectors with components of
%          position and velocity for numDims in {1,2,3}. For example, in
%          3D, the components are [x;y;z;xDot;yDot;zDot].
% useHalfRange A boolean value specifying whether the bistatic range value
%          (and thus the range rate) should be divided by two. This
%          normally comes up when operating in monostatic mode, so that the
%          range reported is a one-way range. The default if this parameter
%          is not provided is false.
%      lTx The (2*numDims)X1 position and velocity vector of the
%          transmitter.
%      lRx The (2*numDims)X1 position and velocity vector of the
%          receiver.
%
%OUTPUTS: JTotal A 1X(2*numDims)XN set of gradients of the bistatic range
%           rates implied by xG with partial derivatives taken with respect
%           to position and velocity components. In 3D, the components are
%           taken in the order [x,y,z,xDot,yDot,zDot]; in 2D and 1D they
%           are taken in the order [x,y,xDot,yDot] and [x,xDot].
%
%The non-relativistic range-rate approximation is derived in [1], where
%expressions for the gradient are given in Appendix I.
%
%EXAMPLE 1:
%Here, we demonstrate that the gradient produced by this function is
%consistent with what one obtains using numerical differentiation.
% lRx=[4e3;-6e3;12;-80;-80;2];
% lTx=[-12e3;8e3;5e3;10;-3;2];
% t=[-3e3;-2e3;-1e3;300;-100;1];
% useHalfRange=false;
% gradVal=rangeRateGradient(t,useHalfRange,lTx,lRx);
% fun=@(t)getRangeRate(t,useHalfRange,lTx,lRx);
% gradNumDiff=numDiff(t,fun,1);
% RelDiff=max(abs((gradNumDiff-gradVal)./gradNumDiff))
%The relative difference between the values should be less than 1e-8, which
%indicates good agreement.
%
%EXAMPLE 2:
%This is similar to example 1, exept the target is the transmitter.
% lRx=[4e3;-6e3;12;-80;-80;2];
% t=[-3e3;-2e3;-1e3;300;-100;1];
% useHalfRange=false;
% gradVal=rangeRateGradient(t,useHalfRange,t,lRx);
% fun=@(t)getRangeRate(t,useHalfRange,t,lRx);
% gradNumDiff=numDiff(t,fun,1);
% RelDiff=max(abs((gradNumDiff-gradVal)./gradNumDiff))
%The relative difference between the values should be around 1e-8, which
%indicates good agreement.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xG,1);
numDim=xDim/2;
N=size(xG,2);

if(nargin<2||isempty(useHalfRange))
   useHalfRange=false; 
end

if(nargin<3||isempty(lTx))
    lTx=zeros(xDim,1); 
end

if(nargin<4||isempty(lRx))
    lRx=zeros(xDim,1); 
end

JTotal=zeros(1,2*numDim,N);
for curPoint=1:N
    x=xG(:,curPoint);
    switch(numDim)
        case 1
            J=zeros(1,2);

            dtr=x(1)-lRx(1);
            dtl=x(1)-lTx(1);

            %The position derivatives.
            J(1)=0;

            %The velocity derivatives.
            J(2)=dtr/norm(dtr);

            if(dtl~=0)
                J(2)=J(2)+dtl/norm(dtl);
            end
        case 2
            J=zeros(1,4);

            dtr=x(1:2)-lRx(1:2);
            dtl=x(1:2)-lTx(1:2);
            dvr=x(3:4)-lRx(3:4);
            dvl=x(3:4)-lTx(3:4);

            drB=dtr/norm(dtr);
            %The position derivatives.
            J(1:2)=A2(dtr)*dvr/norm(dtr)^3;

            if(any(dtl~=0))
                drB=drB+dtl/norm(dtl);
                J(1:2)=vec(J(1:2))+A2(dtl)*dvl/norm(dtl)^3;
            end

            %The velocity derivatives.
            J(3:4)=drB;
        case 3
            J=zeros(1,6);

            dtr=x(1:3)-lRx(1:3);
            dtl=x(1:3)-lTx(1:3);
            dvr=x(4:6)-lRx(4:6);
            dvl=x(4:6)-lTx(4:6);

            drB=dtr/norm(dtr);
            %The position derivatives.
            J(1:3)=A(dtr)*dvr/norm(dtr)^3;

            if(any(dtl~=0))
                drB=drB+dtl/norm(dtl);
                J(1:3)=vec(J(1:3))+A(dtl)*dvl/norm(dtl)^3;
            end

            %The velocity derivatives.
            J(4:6)=drB;
        otherwise
            error('Invalid state dimensionality.')
    end
    JTotal(1,:,curPoint)=J;
end

if(useHalfRange)
    JTotal=JTotal/2; 
end
end

function val=A(vect)
%%A This is a simple helper function as defined in Appendix I of [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.

    x=vect(1);
    y=vect(2);
    z=vect(3);
    val=[y^2+z^2,-x*y,     -x*z;
        -x*y,     x^2+z^2, -y*z;
        -x*z,    -y*z,      x^2+y^2];
end

function val=A2(vect)
%%A2 This is a 2D version of the helper function A as defined in Appendix
%     I of [1].
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.

    x=vect(1);
    y=vect(2);
    val=[y^2, -x*y;
        -x*y,  x^2];
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
