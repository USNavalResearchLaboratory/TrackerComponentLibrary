function [zPred,Pzz,otherInfo]=QMCKalMeasPred(xPred,PPred,zDim,h,numSamples,innovTrans,measAvgFun,stateDiffTrans,stateTrans)
%%QMCKALMEASPRED Perform the measurement prediction part of the measurement
%           update step of the quasi-Monte Carlo Kalman filter with
%           additive measurement noise, as described in Section III of [1],
%           except random samples are used instead of a quasi-random
%           sequence. The function QMCKalUpdateWithPred can be used to
%           complete the measurement update. Separating the measurement
%           prediction step from the rest of the update step can make the
%           creation of multiple measurement association hypotheses from a
%           single target prediction more efficient. The full measurement
%           update function is QMCKalUpdate.
%
%INPUTS: xPred The xDim X 1 predicted target state.
%        PPred The xDim X xDim predicted state covariance matrix.                     
%         zDim The dimensionality of the measurement vector.
%            h A function handle for the measurement function that takes
%              the state as its argument.
%   numSamples The number of samples to use for the Monte Carlo integration
%              in the filter. If this parameter is omitted or an empty
%              matrix is passed, then a default of 100 is used.
%   innovTrans An optional function handle that transforms the value of the
%              difference between the observation and any predicted points.
%              This must be able to handle sets of differences. For a zDim
%              measurement, this must be able to handle a zDimXN matrix of
%              N differences. This only needs to be supplied when a
%              measurement difference must be restricted to a certain
%              range. For example, the innovation between two angles will
%              be 2*pi if one angle is zero and the other 2*pi, even though
%              they are the same direction. In such an instance, a function
%              handle to the wrapRange function with the appropriate
%              parameters should be passed for innovTrans.
%   measAvgFun An optional function handle that, when given N measurement
%              values with weights, produces the weighted average. This
%              function only has to be provided if the domain of the
%              measurement is not linear. For example, when averaging
%              angular values, then the function meanAng should be used.
% stateDiffTrans An optional function handle that, like innovTrans does for
%              the measurements, a difference between states and transforms
%              it however might be necessary. For example, a state
%              containing angular components will generally need to be
%              transformed so that the difference between the angles is
%              wrapped to -pi/pi.
%   stateTrans An optional function that takes a state estimate and
%              transforms it. This is useful if one wishes the elements of
%              the state to be bound to a certain domain. For example, if
%              an element of the state is an angle, one might generally
%              want to bind it to the region +/-pi.
%
%OUTPUTS: zPred The zDimX1 measurement prediction from the filter.
%        PzPred The zDimXzDim covariance matrix associated with zPred.
%     otherInfo A structure containing members of intermediate results of
%               this function that can be passed to QMCKalUpdateWithPred
%               when updating with a measurement.
%
%See the comments to QMCKalUpdate for more information on the algorithm.
%
%EXAMPLE:
%With this example, we demonstrate that one gets the same result using
%QMCKalUpdate in one step as with using QMCKalMeasPred followed by
%QMCKalUpdateWithPred.
% xPred=[1e3;-2e3;100;200];
% h=@(x)([sum(x.^2);x(1)-x(4)^(3/2)]);
% PPred=[28,   3.5,    6,  8.5;
%       3.5,    23,  8.5,   11;
%         6,   8.5,   18, 13.5;
%       8.5,    11, 13.5,   13];
% z=1e6*[5.050000548964568;
%       -0.001829553054023];
% zDim=size(z,1);
% R=eye(zDim,zDim);%Measurement covariance matrix.
% %We use the same random number generator seed for both approaches so that
% %the random part of the results is identical.
% rng('default');
% rng(1);
% %The update in one step.
% [xUpdate,PUpdate,innov,Pzz,W]=QMCKalUpdate(xPred,PPred,z,R,h);
% %Reset the random number generator.
% rng(1);
% %The update in two steps.
% [zPred,PzPred,otherInfo]=QMCKalMeasPred(xPred,PPred,zDim,h);
% [xUpdate1,PUpdate1,innov1,Pzz1,W1]=QMCKalUpdateWithPred(z,R,zPred,PzPred,otherInfo);
% %One will see that the one and two step updates agree.
% max(abs([xUpdate1(:)-xUpdate(:);PUpdate1(:)-PUpdate(:);innov1(:)-innov(:);Pzz1(:)-Pzz(:);W1(:)-W(:)]))
%
%REFERENCES:
%[1] D. Guo and X. Wang, "Quasi-Monte-Carlo filtering in nonlinear
%    dynamical systems," IEEE Transactions on Signal Processing, vol. 54,
%    no. 6, pp. 2087-2098, Jun. 2006.
%
%June 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(numSamples))
    numSamples=100; 
end

if(nargin<6||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(x)x;
end

if(nargin<7||isempty(measAvgFun))
    measAvgFun=@(zPoints,w)calcMixtureMoments(zPoints,w);
end

if(nargin<8||isempty(stateDiffTrans))
   stateDiffTrans=@(x)x; 
end

if(nargin<9||isempty(stateTrans))
    stateTrans=@(x)x; 
end

xDim=size(xPred,1);

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

%Centered, predicted cubature measurement points, transformed as
%necessary to keep the values within a desired range.
zPredCenPoints=innovTrans(bsxfun(@minus,zPredPoints,zPred));

%Here, R is included in Pzz, rather than having to be added separately as
%in Equation 19.
Pzz=zeros(zDim,zDim);
Pxz=zeros(xDim,zDim);
for curP=1:numSamples
    diff=zPredCenPoints(:,curP);
    Pzz=Pzz+w*(diff*diff');
    Pxz=Pxz+w*stateDiffTrans(xSamples(:,curP)-xPred)*diff';
end
%Pxz is not needed for the measurement prediction, but we compute it
%here, so that it need not be recomputed again and again if
%QMCKalUpdateWithPred is called for multiple measurements.

otherInfo.innovTrans=innovTrans;
otherInfo.stateTrans=stateTrans;
otherInfo.xPred=xPred;
otherInfo.PPred=PPred;
otherInfo.Pxz=Pxz;

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
