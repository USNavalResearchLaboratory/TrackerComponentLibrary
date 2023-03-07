function [mu,P]=calcMixtureMoments(xi,w,PHyp,muHyp,diffTransFun,numMergeDims,omitSpreadTerms)
%%CALCMIXTUREMOMENTS Given a set of sample vectors xi and associated
%                    weights w, calculate the first two moments of the
%                    mixture. Those are, the mean mu and the covariance
%                    matrix P. An additional input can allow the function
%                    to be used for computing the covariance matrix of a
%                    Gaussian mixture, when xi is the set of mean vectors,
%                    w is the set of  weights and the third input PHyp is
%                    the set of covariance matrices for the hypotheses.
%                    The input muHyp makes P on the output the mean squared
%                    error matrix for an estimate muHyp rather than a
%                    covariance matrix.
%
%INPUTS: xi A numDimXnumPoints matrix of column vectors.
%         w A numPointsX1 or 1XnumPoints vector of the weights associated
%           with each of the column vectors. Normally all w>0 to assure
%           that P is positive definite. It is assumed that the weights sum
%           to 1 when computing the mu and P. If this parameter is omitted
%           or an empty matrix is passed, then all values in w are set to
%           1/numPoints (uniform distribution).
%      PHyp An optional numDimXnumDimXnumPoints parameter. If one wants to
%           find the covariance matrix of a mixture distribution, then this
%           is the set of covariances of the individual components (whose
%           means are xi). If omitted or an empty matrix is passed, P will
%           be computed as the covariance matrix of a bunch of points. If P
%           is not requested on the output, this argument has no effect.
%     muHyp An optional numDimX1 input that changes how the covariance
%           matrix is computed. If something other than an empty matrix is
%           provided, the output P will not be computed centered about mu
%           (it will not be the covariance) but will be computed centered
%           around muHyp, thus making P a mean-squared error matrix. This
%           input can be useful when using calcMixtureMoments to compute
%           updates for the global nearest neighbor joint probabilistic
%           data association filter.
% diffTransFun An optional function handle that transforms differences
%           between values when computing the covariance. Though
%           calcMixtureMoments is not meant for use with circular data,
%           passing a function to wrap data can be used as an ad-hoc fix
%           for computing a covariance matrix when using circular data,
%           such as when dealing with longitudinal values.
% numMergeDims It is possible that numDim> the number of dimensions that
%           one actually cares about for the merged value mu. If so, then
%           this is the number of dimensions that should actually be
%           merged. This discrepancy can arise when reducing models in an
%           IMM for display where some models have additional non-mixing
%           components. The default if this parameter is omitted or an
%           empty matrix is passed is numDim.
% omitSpreadTerms An optional term that if ture makes P the average of PHyp
%           rather than the true covariance matrix of the mixture. The
%           default if omitted or an empty matrix is passed is false. As in
%           Chapter 1.4.16 of [1], the covariance matrix computation
%           includes a spread of the means term. If this is true, then that
%           term is omitted.
%
%OUTPUTS: mu The numMergeDimsX1 mean of the mixture.
%          P The numMergeDimsXnumMergeDims covariance matrix of the
%            mixture.
%
%This function can be used to conveniently compute the mean and covariance
%of a set of cubature points and weights. It can also be used to compute
%the mean and covariance of a Gaussian mixture distribution, if xi are the
%means of the components and w are the weights associated with each
%component in the mixture. This function is also useful in reducing
%mixtures in multiple model algorithms, and reducing hypotheses in data
%association algorithms.
%
%Equations for computing the moments of mixtures are derived in Chapter
%1.4.16 of [1].
%
%The most common uses of this function are of the form
%[mu,P]=calcMixtureMoments(xi,w) for computing the mean and covariance
%matrix of a batch of weighted samples and also of the form 
%[mu,P]=calcMixtureMoments(xi,w,PHyp) for computing the mean and covariance
%matrix of a Gaussian mixture distribution.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kirubarajan, Estimation with
%    Applications to Tracking and Navigation. New York: John Wiley and
%    Sons, Inc, 2001.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(omitSpreadTerms))
    omitSpreadTerms=false;
end

if(nargin<6||isempty(numMergeDims))
    numMergeDims=size(xi,1);
end

if(nargin<5)
    diffTransFun=[];
end

if(nargin<4)
    muHyp=[];
end

if(nargin<3)
    PHyp=[];
end

numPoints=size(xi,2);

if(nargin<2||isempty(w))
    w=ones(numPoints,1);
end

%Make sure that w is a normalized column vector
w=w(:)/sum(w(:));

mu=sum(bsxfun(@times,xi(1:numMergeDims,:),w.'),2);%The mean

if(nargout>1)
    if(~isempty(muHyp))
        mean2Use=muHyp(1:numMergeDims);
    else
        mean2Use=mu;
    end
    
    P=zeros(numMergeDims,numMergeDims);

    if(omitSpreadTerms==false)
        if(~isempty(diffTransFun))
            diff=diffTransFun(bsxfun(@minus,xi(1:numMergeDims,:),mean2Use));
        else
            diff=bsxfun(@minus,xi(1:numMergeDims,:),mean2Use);
        end
    
        if(~isempty(PHyp))%If covariance matrices are provided.
            for curPoint=1:numPoints
                P=P+w(curPoint)*(PHyp(1:numMergeDims,1:numMergeDims,curPoint)+diff(:,curPoint)*diff(:,curPoint)');
            end
        else
            for curPoint=1:numPoints
                P=P+w(curPoint)*(diff(:,curPoint)*diff(:,curPoint)');
            end
        end
    elseif(~isempty(PHyp))
        %If covariance matrices are provided and we are not including the
        %spread terms (A solution that is generally incorrect).
        for curPoint=1:numPoints
            P=P+w(curPoint)*PHyp(1:numMergeDims,1:numMergeDims,curPoint);
        end
    end
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
