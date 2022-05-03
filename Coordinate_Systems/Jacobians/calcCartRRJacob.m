function J=calcCartRRJacob(components,xState,useHalfRange,lTx,lRx)
%%CALCCARTRRJACOB Compute the Jacobian matrix for a Cartesian measurement
%           possibly with a bistatic range rate component in 1D, 2D or 3D
%           for a Cartesian target state. This function can be useful when
%           using an extended Kalman filter (EKF) with Cartesian-converted
%           measurements and range rate. Without the range rate component,
%           the result is just the measurement matrix, which, when
%           multiplied by the target state, extracts the position
%           components and is often designated by H in the discrete-time
%           Kalman filter.
%
%INPUTS: components This specified whether the output should contain a
%               range-rate component. Possible values are:
%               0 The output contains position and range rate rows. This is
%                 the default if an empty matrix is passed.
%               1 The output just has position rows.
%        xState The xDimX1 target state vector in the global coordinate
%               system with [x;y;z;xDot;yDot;zDot] components in 3D or
%               [x;y;xDot;yDot] components in 2D.
%  useHalfRange A boolean value specifying whether the bistatic range
%               (and thus the range rate) value has been divided by two.
%               This normally comes up when operating in monostatic mode,
%           lTx The transmitter state vector in the global coordinate
%               system with [x;y;z;xDot;yDot;zDot] components in 3D or 
%               [x;y;xDot;yDot] components in 2D. If omitted, then a vector
%               of zeros is used.
%           lRx The receiver state vector in the global coordinate system
%               with [x;y;z;xDot;yDot;zDot] components in 3D or
%               [x;y;xDot;yDot] components in 2D. If omitted, then a vector
%               of zeros is used.
%
%OUTPUTS: J The (xDim/2+1)XxDim Jacobian matrix. Each row is a component
%           of x-y-z-range rate and each column is the derivative with
%           respect to x,y,z,xDot,yDot,zDot. If range rate is not desired,
%           then the matrix is (xDim/2)XxDim in size.
%
%The derivatives of the position components with themselves are clear (1
%for common components and zero for different components). The derivatives
%of the range rate components are obtained using the rangeRateGradient
%function.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(xState,1);
posDim=xDim/2;
measDim=posDim+1;

if(isempty(components))
    components=0;
end

if(components==0)
    if(nargin<3||isempty(useHalfRange))
        useHalfRange=false;
    end

    if(nargin<4||isempty(lTx))
        lTx=zeros(xDim,1); 
    end

    if(nargin<5||isempty(lRx))
        lRx=zeros(xDim,1); 
    end

    J=zeros(measDim,xDim);
    %Position components are just passed through.
    J(1:posDim,1:posDim)=eye(posDim,posDim);

    J(posDim+1,:)=rangeRateGradient(xState,useHalfRange,lTx,lRx);
else
    J=zeros(posDim,xDim);
    J(1:posDim,1:posDim)=eye(posDim,posDim);
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
