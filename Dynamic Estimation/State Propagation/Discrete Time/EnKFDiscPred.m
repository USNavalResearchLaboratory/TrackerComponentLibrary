function [xEnsemb,xPred,PPred,vSamp]=EnKFDiscPred(xEnsemb,f,SQ,filterType,stateDiffTrans,vSamp)
%%ENKFDISCPRED Perform the state prediction step of a discrete-time
%           ensemble Kalman filter. Such a filter is similar to the pure
%           propagation filter, but adjusts the covariance of the
%           propagated points by adding random noise rather than scaling or
%           adding deterministic values to achieve a desired magnitude
%           increase in the sample covariance. The input and primary output
%           of this function is a set of unweighted sample points. Thus,
%           information on higher order moments is, to an extent,
%           preserved. The filter propagates the points, which can then be
%           given to the measurement update function. The first two moments
%           (mean and covariance) are also available for display, but do
%           not directly play a role in the filter.
%
%INPUTS: xEnsemb  An xDim X numSamples matrix of the target state samples
%                 (ensemble points) that are to be predicted.
%          f A function handle for the state transition function that takes
%            the state as its parameter. If the process noise is not
%            additive, then it takes the state and the process noise as
%            parameters.
%         SQ The xDimX xDim lower-triangular square root of the process
%            noise covariance matrix. If vSamp is provided, then an empty
%            matrix can be passed for this.
% filterType A parameter effecting how the filter performs the prediction.
%            Possible values are:
%            0) (The default if omitted or an empty matrix is passed) The
%               process noise is additive.
%            1) The process noise is not additive.
% stateDiffTrans An optional function handle that takes an xDimXN matrix
%                of N differences between states estimates and
%                transforms them however might be necessary. For
%                example, a state continaing angular components will
%                generally need differences between angular components
%                wrapped to the range +/-pi. This is only needed if PPred
%                is desired to be returned.
%      vSamp An optional xDimXnumSamples matrix of noise samples that are
%            to be used to perturb the predicted states due to the process
%            noise. This is available as a parameter in case one wishes to
%            use a deterministic sampling method. If omitted, the noise
%            samples are randomly drawn from a N(O,SQ*SQ') Gaussian
%            distribution and forced to be zero-mean. This is consistent
%            with the suggestion in [2] that strives to avoid changing the
%            best estimates.
%
%The implementation is consistent with the descriptions given in [1] and
%[2].
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

if(nargin<4||isempty(filterType))
    filterType=0;
end

if(nargin<5||isempty(stateDiffTrans))
    stateDiffTrans=@(x)x;
end

xDim=size(xEnsemb,1);
numSamples=size(xEnsemb,2);

%Generate the noise samples, if no samples were given.
if(nargin<6||isempty(vSamp))
    vSamp=SQ*randn(xDim,numSamples);
    
    %Consistent with [2], force the random samples to be zero-mean.
    vMean=mean(vSamp,2);
    vSamp=bsxfun(@minus,vSamp,vMean);
end

%Propagate the state samples
if(filterType==0)
    %Additive noise.
    for curP=1:numSamples
        xEnsemb(:,curP)=f(xEnsemb(:,curP))+vSamp(:,curP);
    end
elseif(filterType==1)
    %Non-additive noise.
    for curP=1:numSamples
        xEnsemb(:,curP)=f(xEnsemb(:,curP),vSamp(:,curP));
    end
else
end

if(nargout>1)
    xPred=stateAvgFun(xEnsemb,2);
    xPredCenPoints=stateDiffTrans(bsxfun(@minus,xEnsemb,xPred));%Center the samples.
    %The updated covariance
    PPred=(1/(numSamples-1))*(xPredCenPoints*xPredCenPoints');
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
