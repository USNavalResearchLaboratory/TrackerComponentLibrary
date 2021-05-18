function [xUpdate,PUpdate,innov,Pzz,W]=QMCKalUpdate(xPred,PPred,z,R,h,numSamples,innovTrans,measAvgFun,stateDiffTrans,stateTrans)
%%QMCKALUPDATE Perform the measurement update step in the quasi-Monte Carlo
%              Kalman filter with additive measurement noise, as described
%              in Section III of [1], except random samples are used
%              instead of a quasi-random sequence. Note that the non-
%              quadratic form of the covariance update does not 100%
%              guarantee that the resulting covariance matrix will be
%              positive (semi-)definite. This filter is essentially the
%              cubature Kalman filter update using random cubature points.
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        PPred The xDim X xDim predicted state covariance matrix.                     
%            z The zDim X 1 vector measurement.
%            R The zDim X zDim measurement covariance matrix in the native
%              coordinate system of the measurement.
%            h A function handle for the measurement function that takes
%              the state as its argument.
%   numSamples The number of samples to use for the Monte Carlo integration
%              in the filter. If this parameter is omitted or an empty
%              matrix is passed, then a default of 100 is used.
%   innovTrans An optional function handle that computes and optionally
%              transforms the value of the difference between the
%              observation and any predicted points. This is called as
%              innovTrans(a,b) and the default if omitted or an empty
%              matrix is passed is @(a,b)bsxfun(@minus,a,b). This must be
%              able to handle sets of values. For a zDimX1 measurement,
%              either of the inputs could be zDimXN in size while one of
%              the inputs could be zDimX1 in size.  This only needs to be
%              supplied when a measurement difference must be restricted
%              to a certain range. For example, the innovation between two
%              angles will be 2*pi if one angle is zero and the other
%              2*pi, even though they are the same direction. In such an
%              instance, a function handle to the
%              wrapRange(bsxfun(@minus,a,b),-pi,pi) function with the
%              appropriate parameters should be passed for innovTrans.
%   measAvgFun An optional function handle that, when given N measurement
%              values with weights, produces the weighted average. This
%              function only has to be provided if the domain of the
%              measurement is not linear. For example, when averaging
%              angular values, then the function meanAng should be used.
% stateDiffTrans An optional function handle that takes an xDimXN matrix of
%              N differences between states and transforms them however
%              might be necessary. For example, a state containing angular
%              components will generally need to be transformed so that the
%              difference between the angles is wrapped to -pi/pi.
%   stateTrans An optional function that takes a state estimate and
%              transforms it. This is useful if one wishes the elements of
%              the state to be bound to a certain domain. For example, if
%              an element of the state is an angle, one might generally
%              want to bind it to the region +/-pi.
%
%OUTPUTS: xUpdate The xDim X 1 updated state vector.
%         PUpdate The updated xDim X xDim state covariance matrix.
%      innov, Pzz The zDimX1 innovation and the zDimXzDim innovation
%                 covariance matrix are returned in case one wishes to
%                 analyze the consistency of the estimator or use those
%                 values in gating or likelihood evaluation.
%               W The xDimXzDim gain used in the update. This can be
%                 useful when gating and using the function
%                 calcMissedGateCov.
%
%This implements the measurement update step in Section III of [1]
%(Algorithm 1), but here we use random numbers instead of a quasi-random
%sequence. This is essentially a Cubature Kalman filter that performs
%Monte Carlo integration rather than cubature integration.
%
%The optional parameters innovTrans and measAvgFun are not described in
%references [1], but allow for possible modifications to the filter as
%described in [2]. The parameters have been added to allow the filter to be
%used with angular quantities. For example, if the measurement consisted
%of range and angle, z=[r;theta], then
%innovTrans=@(a,b)[bsxfun(@minus,a(1,:),b(1,:));
%                  wrapRange(bsxfun(@minus,a(2,:),b(2,:)),-pi,pi)];
%measAvgFun=@(z,w)[calcMixtureMoments(z(1,:),w);
%                  meanAng(z(2,:),w')];
%should be used to approximately deal with the circular nature of the
%measurements.
%
%REFERENCES:
%[1] D. Guo and X. Wang, "Quasi-Monte-Carlo filtering in nonlinear
%    dynamical systems," IEEE Transactions on Signal Processing, vol. 54,
%    no. 6, pp. 2087-2098, Jun. 2006.
%[2] D. F. Crouse, "Cubature/ unscented/ sigma point Kalman filtering with
%    angular measurement models," in Proceedings of the 18th International
%    Conference on Information Fusion, Washington, D.C., 6-9 Jul. 2015.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(numSamples))
   numSamples=100; 
end

if(nargin<7||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(a,b)bsxfun(@minus,a,b);
end

if(nargin<8||isempty(measAvgFun))
   measAvgFun=@(zPoints,w)calcMixtureMoments(zPoints,w);
end

if(nargin<9||isempty(stateDiffTrans))
   stateDiffTrans=@(x)x; 
end

if(nargin<10||isempty(stateTrans))
   stateTrans=@(x)x; 
end

xDim=size(xPred,1);
zDim=size(z,1);

%Generate the samples from the prediction.
xSamples=bsxfun(@plus,xPred,cholSemiDef(PPred,'lower',1)*randn(xDim,numSamples));
w=1/numSamples;%uniform weighting.

%Predicted measurement points
zPredPoints=zeros(zDim,numSamples);
for curP=1:numSamples
    zPredPoints(:,curP)=h(xSamples(:,curP));
end

%Measurement prediction.
zPred=measAvgFun(zPredPoints,repmat(1/numSamples,[numSamples,1]));

%The innovation, transformed as necessary to keep values in a desired
%range.
innov=innovTrans(z,zPred);

%Here, R is included in Pzz, rather than having to be added separately as
%in Equation 19.
Pzz=R;
Pxz=zeros(xDim,zDim);
for curP=1:numSamples
    diff=innovTrans(zPredPoints(:,curP),zPred);
    Pzz=Pzz+w*(diff*diff');
    Pxz=Pxz+w*stateDiffTrans(xSamples(:,curP)-xPred)*diff';
end

%The filter gain
W=Pxz/Pzz;

%Updated state estimate
xUpdate=stateTrans(xPred+W*innov);

%Updated state covariance matrix
PUpdate=PPred-W*Pzz*W';
%Ensure symmetry
PUpdate=(PUpdate+PUpdate')/2;
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
