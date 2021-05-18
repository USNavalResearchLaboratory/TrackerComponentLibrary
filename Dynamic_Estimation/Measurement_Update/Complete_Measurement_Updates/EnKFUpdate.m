function [xEnsemb,xUpdate,PUpdate,wSamp,Pzz,W,innovPoints]=EnKFUpdate(xEnsemb,z,SR,h,filterType,innovTrans,measAvgFun,stateDiffTrans,stateAvgFun,stateTrans,wSamp)
%ENKFUPDATE Perform the measurement update step in the ensemble Kalman
%           filter. Such a filter is similar to the pure propagation
%           filter, but adjusts the covariance of the propagated points by
%           adding random noise rather than scaling or adding deterministic
%           values to achieve a desired magnitude increase in the sample
%           covariance. The input and primary output of this function is a
%           set of unweighted sample points. Thus, information on higher
%           order moments is, to an extent, preserved. The filter updates
%           the points, which can then be given to the propagation
%           function. The first two moments (mean and covariance) are also
%           available for display, but do not directly play a role in the
%           filter.
%
%INPUTS: xEnsemb An xDim X numSamples matrix of the prior predicted target
%                state samples (ensemble points).
%              z The zDim X 1 vector measurement.
%             SR The zDim X zDim lower-triangular square root of the
%                measurement covariance matrix. If wSamp is provided, then
%                an empty matrix can be passed for this.
%              h A function handle for the measurement function that takes
%                the state as its argument (z=h(x)). If filterType=1, then
%                h takes the state and the measurement noise as its
%                arguments (z=h(x,w))
%     filterType A parameter effecting how the filter performs the
%                measurement update. Possible values are:
%                0) (The default if omitted or an empty matrix is passed)
%                   The randomly generated (additive) measurement noise is
%                   added to the predicted samples prior to computing the
%                   cross correlations for the gain. This is not consistent
%                   with [1] and [2], but is consistent with the normal
%                   definition of the innovation covariance used in the
%                   gain.
%                1) This is the same as filterType 0, but with non-
%                   additive noise, so h has the form h(x,w), where w is
%                   the measurement noise.
%                2) Here, the filter is implemented as in [1] and [2],
%                   where the additive noise is not included in the
%                   computation of the gain. Rather, it is only included in
%                   the innovation term.
%    innovTrans An optional function handle that computes and optionally
%               transforms the value of the difference between the
%               observation and any predicted points. This is called as
%               innovTrans(a,b) and the default if omitted or an empty
%               matrix is passed is @(a,b)bsxfun(@minus,a,b). This must be
%               able to handle sets of values. For a zDimX1 measurement,
%               either of the inputs could be zDimXN in size while one of
%               the inputs could be zDimX1 in size.  This only needs to be
%               supplied when a measurement difference must be restricted
%               to a certain range. For example, the innovation between two
%               angles will be 2*pi if one angle is zero and the other
%               2*pi, even though they are the same direction. In such an
%               instance, a function handle to the
%               wrapRange(bsxfun(@minus,a,b),-pi,pi) function with the
%               appropriate parameters should be passed for innovTrans.
%     measAvgFun An optional function handle that, when given N measurement
%                values, produces the average. This function only has to be
%                provided if the domain of the measurement is not linear.
%                For example, when averaging angular values, then the
%                function meanAng should be used.
% stateDiffTrans An optional function handle that, like innovTrans does for
%                the measurements, takes an xDimXN matrix of N differences
%                between states and transforms them however might be
%                necessary. For example, a state containing angular
%                components will generally need to be transformed so that
%                the difference between the angles is wrapped to -pi/pi.
%    stateAvgFun An optional function that given an xDimXN matrix of N
%                state estimates and provides the average of the state
%                estimates. This is necessary if, for example, states with
%                angular components are averaged.
%     stateTrans An optional function that takes a matrix of N state
%                estimates and transforms them. This is useful if one
%                wishes the elements of the state to be bound to a certain
%                domain. For example, if an element of the state is an
%                angle, one might generally want to bind it to the region
%                +/-pi.
%          wSamp An optional zDimXnumSamples matrix of noise samples that
%                are to be used to perturb measurement terms. This is
%                available as a parameter in case one wishes to use a
%                deterministic sampling method. If omitted, the noise
%                samples are randomly drawn from an N(O,SR*SR') Gaussian
%                distribution and forced to be zero-mean. This is
%                consistent with the suggestion in [2] that strives to
%                avoid changing the best estimates.
%
%OUTPUTS: xEnsemb The updated ensemble points.
%         xUpdate The mean of the updated ensemble points.
%         PUpdate The covariance matrix associated with the updated sample
%                 points.
%           wSamp The measurement noise samples used. If wSamp is provided
%                 as an input, this is the input. Otherwise, this is the
%                 randomly generated samples (with forced zero-mean).
% Pzz, innovPoints The zDimXzDim innovation covariance matrix and the
%                 zDimXnumSamples innovation value for each of the points
%                 are returned in case one wishes to analyze the
%                 consistency of the estimator or use those values in
%                 gating or likelihood evaluation. One will typically use
%                 the average of innovPoints over all of the samples.
%               W The xDimXzDim gain used in the update. This can be
%                 useful when gating and using the function
%                 calcMissedGateCov.
%
%The basic algorithm is described in [1] and [2]. However, as described,
%the innovation covariance (the covariance matrix of the measurement and
%the prediction) would be underestimated as the matrix is computed without
%the measurement noise added. Thus, here we offer multiple implementations.
%For filterType=0, the measurement noise terms are added to the predicted
%measurements prior to computing the innovation covariance, which woul
%would assume is the correct way. Such a method can be easily modified to
%handle non-additive noise giving us filterType=1. On the other hand, in
%filterType=2, the terms are added only to the measurement as in [1] and
%[2].
%
%Consistent with the suggestion in [2], rather than adding truely
%Gaussian noise, the randomly chosen sample is modified to force it to have
%zero mean (if wSamp is not provided).
%
%For additive noise, the measurement equation is z=h(x)+w, where w is
%noise. For non-additive noise, it is z=f(x,w).
%
%REFERENCES:
%[1] S. Gillijns, O. B. Mendoza, J. Chandrasekar, B. L. R. De Moor, D. S.
%    Bernstein, and A. Ridley, "What is the ensemble Kalman filter and how
%    well does it work?" in Proceedings of the 2006 American Control
%    Conference, Minneapolis, MN, 14-16 Jun. 2006, pp. 4448-4453.
%[2] P. L. Houtekamer and H. L. Mitchell, "Ensemble Kalman filtering,"
%    Quarterly Journal of the Royal Meteorological Society, vol. 131, no.
%    613, pp. 3269-3289, Oct. 2005.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

zDim=size(z,1);
numSamples=size(xEnsemb,2);

if(nargin<5||isempty(filterType))
   filterType=0; 
end

if(nargin<6||isempty(innovTrans))
    %The function just returns the input.
    innovTrans=@(a,b)bsxfun(@minus,a,b);
end

if(nargin<7||isempty(measAvgFun))
    measAvgFun=@(zPoints)mean(zPoints,2);
end

if(nargin<8||isempty(stateDiffTrans))
    stateDiffTrans=@(x)x; 
end

if(nargin<9||isempty(stateAvgFun))
    stateAvgFun=@(xPoints)mean(xPoints,2);
end

if(nargin<10||isempty(stateTrans))
    stateTrans=@(x)x; 
end

%Generate the noise samples, if no samples were given.
if(nargin<11||isempty(wSamp))
    wSamp=SR*randn(zDim,numSamples);
    
    %Consistent with [2], force the random samples to be zero-mean.
    wMean=mean(wSamp,2);
    wSamp=bsxfun(@minus,wSamp,wMean);
end

%Allocate space for the perturbed observations
zPertPoints=zeros(zDim,numSamples);

%Get the perturbed measurements based on the filter type chosen.
if(filterType==0)
    %Implement in the manner that appears to be correct: the additive noise
    %is added to the predicted values and the innovation is thus taken with
    %respect to the measurement and the perturbed predicted values.
    for curP=1:numSamples
        zPertPoints(:,curP)=h(xEnsemb(:,curP))+wSamp(:,curP);
    end
elseif(filterType==1)
    %Implement in the same manner as filterType 0, but for non-additive
    %noise.
    for curP=1:numSamples
        zPertPoints(:,curP)=h(xEnsemb(:,curP),wSamp(:,curP));
    end
elseif(filterType==2)
    %Implement in the manner consistent with the description in [1] and
    %[2] for additive noise. In this instance, no noise is added to the
    %predicted measurement terms and noise is only added directly to the
    %measurements affecting the innovation in the update.
    for curP=1:numSamples
        zPertPoints(:,curP)=h(xEnsemb(:,curP));
    end
else
     error('Invalid filter type given')
end

zPred=measAvgFun(zPertPoints);

%Centered, perturbed sample measurement points, transformed as
%necessary to keep the values within a desired range.
zPertCenPoints=innovTrans(zPertPoints,zPred);

%The innovation covariance
Pzz=(1/(numSamples-1))*(zPertCenPoints*zPertCenPoints');

xPred=stateAvgFun(xEnsemb);%The mean state
xPredCenPoints=stateDiffTrans(bsxfun(@minus,xEnsemb,xPred));%Center the samples.

%The cross covariance
Pxz=(1/(numSamples-1))*(xPredCenPoints*zPertCenPoints');

%The gain (pinv is used in case the random sampling makes Pzz nearly
%singular).
W=Pxz*pinv(Pzz);

innovPoints=zeros(zDim,1);
if(filterType==0||filterType==1)
    %Update the ensemble points, without noise added to the measurement,
    %because it is included in the predicted points.
    for curP=1:numSamples
        innovPoints(:,curP)=innovTrans(z,zPertPoints(:,curP));
        xEnsemb(:,curP)=stateTrans(xEnsemb(:,curP)+W*innovPoints(:,curP));
    end
else
    %Update the ensemble points adding in the noise to the measurement,
    %because it is not included in the ensemble points. This is how it is
    %in [1] and [2].
    for curP=1:numSamples
        innovPoints(:,curP)=innovTrans(z+wSamp(:,curP),zPertPoints(:,curP));
        xEnsemb(:,curP)=stateTrans(xEnsemb(:,curP)+W*innovPoints(:,curP));
    end
end

if(nargout>1)
    xUpdate=stateAvgFun(xEnsemb);
    xUpdateCenPoints=stateDiffTrans(bsxfun(@minus,xEnsemb,xUpdate));%Center the samples.
    %The updated covariance
    PUpdate=(1/(numSamples-1))*(xUpdateCenPoints*xUpdateCenPoints');
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
