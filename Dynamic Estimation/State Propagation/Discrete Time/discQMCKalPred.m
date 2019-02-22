function [xPred,PPred]=discQMCKalPred(xPrev,PPrev,f,Q,numSamples,stateDiffTrans,stateAvgFun,stateTrans)
%%DISCQMCKALPRED Perform the discrete-time prediction step that comes with
%                the quasi-Monte Carlo Kalman filter with additive process
%                noise.
%
%INPUTS: xPrev The xDim X 1 state estimate at the previous time-step.
%        PPrev The xDim X xDim state covariance matrix at the previous
%              time-step.
%            f A function handle for the state transition function
%              that takes the state as its parameter.
%            Q The xDimX xDim process noise covariance matrix.
%   numSamples The number of samples to use for the Monte Carlo integration
%              in the filter. If this parameter is omitted or an empty
%              matrix is passed, then a default of 100 is used.
% stateDiffTrans An optional function handle that takes an xDimXN matrix of
%              N differences between states estimates and transforms them
%              however might be necessary. For example, a state continaing
%              angular components will generally need differences between
%              angular components wrapped to the range +/-pi.
%  stateAvgFun An optional function that given an xDimXN matrix of N state
%              estimates and an NX1 vector of weights, provides the
%              weighted average of the state estimates. This is necessary
%              if, for example, states with angular components are
%              averaged.
%   stateTrans An optional function that takes a state estimate and
%              transforms it. This is useful if one wishes the elements of
%              the state to be bound to a certain domain. For example, if
%              an element of the state is an angle, one should generally
%              want to bind it to the region +/-pi. This is not applied to
%              the output of f.
%
%OUTPUTS: xPred The xDimX1 predicted state estimate.
%         PPred The xDimXxDim predicted state covariance estimate.
%
%This implements the prediction step in Section III of [1] (Algorithm 1).
%This is essentially a Cubature Kalman filter that performs Monte Carlo
%integration rather than cubature integration.
%
%REFERENCES:
%[1] D. Guo and X. Wang, "Quasi-Monte-Carlo filtering in nonlinear
%    dynamical systems," IEEE Transactions on Signal Processing, vol. 54,
%    no. 6, pp. 2087-2098, Jun. 2006.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(numSamples))
   numSamples=100; 
end

if(nargin<6||isempty(stateDiffTrans))
    stateDiffTrans=@(x)x;
end

if(nargin<7||isempty(stateAvgFun))
    stateAvgFun=@(x,w)calcMixtureMoments(x,w);
end

if(nargin<8||isempty(stateTrans))
    stateTrans=@(x)x;
end

xDim=size(xPrev,1);

%Generate the samples from the prior.
xSamples=stateTrans(bsxfun(@plus,xPrev,cholSemiDef(PPrev,'lower')*randn(xDim,numSamples)));

%Propagate the samples
for curSamp=1:numSamples
    xSamples(:,curSamp)=f(xSamples(:,curSamp));
end

%Equation 13, calculate the predicted state
xPred=stateAvgFun(xSamples,repmat(1/numSamples,[numSamples,1]));

%Equation 14, calculate the predicted covariance matrix.
xSamples=stateDiffTrans(bsxfun(@minus,xSamples,xPred));
PPred=Q+(1/numSamples)*(xSamples*xSamples');
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
