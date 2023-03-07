function J=rangeGradient(x,useHalfRange,lTx,lRx)
%%RANGEGRADIENT Determine the gradient of a bistatic range measurement (in
%           1D, 2D, 3D) with respect to Cartesian position (gradient
%           components for velocity etc. are not provided). Atmospheric and
%           other propagation effects are not taken into account. The
%           transmitter can be collocated with the target (in which case it
%           is assumed they move together).
%
%INPUTS: x A numPosDimXN set of N target position vectors of the form
%          [x], [x;y] or [x;y;z].
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up
%          when operating in monostatic mode, so that the range reported is
%          a one-way range. The default if this parameter is not provided
%          is false.
%      lTx The numPosDimX1 position vector of the transmitter. If this
%          parameter is omitted or an empty matrix is passed, then the
%          transmitter is assumed to be at the origin.
%      lRx The numPosDimX1 position vector of the receiver. If this
%          parameter is omitted or an empty matrix is passed, then the
%          receiver is assumed to be at the origin.
%
%OUTPUTS: J A 1XnumPosDimXN set of partial derivatives of the bistatic
%           range with respect to components [x,y,z] in 3D, [x,y] in 2D, or
%           [x] in 1D  in that order for each point in x.
%
%Gradients with respect to bistatic range are discussed in [1]. If the
%transmitter is collocated with the target, it is assumed to move with the
%target.
%
%EXAMPLE 1:
%Here, we demonstrate that the gradient produced by this function is
%consistent with what one obtains using numerical differentiation.
% lRx=[4e3;-6e3;12];
% lTx=[-12e3;8e3;5e3];
% t=[-3e3;-2e3;-1e3];
% useHalfRange=false;
% gradVal=rangeGradient(t,useHalfRange,lTx,lRx);
% fun=@(t)getRange(t,useHalfRange,lTx,lRx);
% gradNumDiff=numDiff(t,fun,1);
% RelDiff=max(abs((gradNumDiff-gradVal)./gradNumDiff))
%The relative difference between the values should be less than 1e-9, which
%indicates good agreement.
%
%EXAMPLE 2:
%In this example, we also verify that the results agree with finite
%differencing, but here the transmitter is collocated with the target (and
%thus is assumed to move with the target).
% lRx=[4e3;-6e3;12];
% t=[-3e3;-2e3;-1e3];
% useHalfRange=false;
% gradVal=rangeGradient(t,useHalfRange,t,lRx);
% fun=@(t)getRange(t,useHalfRange,t,lRx);
% gradNumDiff=numDiff(t,fun,1);
% RelDiff=max(abs((gradNumDiff-gradVal)./gradNumDiff))
%The relative difference between the values should again be less than 1e-9,
%which indicates good agreement.
%
%REFERENCES:
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2017 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(x,1);
N=size(x,2);

if(nargin<2||isempty(useHalfRange))
    useHalfRange=false; 
end

if(nargin<4||isempty(lRx))
    lRx=zeros(numDim,1); 
end

if(nargin<3||isempty(lTx))
    lTx=zeros(numDim,1); 
end

J=zeros(1,numDim,N);
for curX=1:N
    deltaRx=x(:,curX)-lRx;
    deltaTx=x(:,curX)-lTx;

    if(all(deltaTx==0))
        %If the transmitter is collocated with the target.
        J(1,:,curX)=deltaRx.'/norm(deltaRx);
    else
        J(1,:,curX)=deltaRx.'/norm(deltaRx)+deltaTx.'/norm(deltaTx);
    end
end

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
