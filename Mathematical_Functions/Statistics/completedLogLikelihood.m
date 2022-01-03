function [CLL,tau] = completedLogLikelihood(data,likelihood,clusterWeights,params,tau)
%%COMPLETEDLOGLIKELIHOOD Computes the conditional expected value of the
%                        completed log-likelihood (CLL) of a cluster model.
%                        This is part of the E step of the expectation
%                        maximization algorithm described in [1]. To
%                        compute the log-likelihood, give tau as an empty
%                        matrix.
%
%INPUTS:
% data: A numDim-by-numSamples matrix where numDim is the state dimension
%       and numSamples is the number of vectors representing individual
%       data points.
% likelihood: A function handle which takes as inputs a single data vector
%             from the input data and a single row of parameters from the
%             input params in that order.
% clusterWeights: A vector of normalizable weights with length equal to the
%                 number of clusters in the model (numClusters).
% params: A numClusters-by-numParams cell array where each row corresponds
%         to the model parameters for the corresponding cluster and
%         params{clust,:} should unpack the model parameters in the correct
%         order to be given to likelihood.
% tau: A numClusters-by-numParams matrix where the i,j entry indicates
%      the jth sample originated from the ith cluster with probability
%      tau(i,j). If left empty, an all ones matrix is used which will cause
%      the returned CLL value to be the traditional log-likelihood. If
%      given but tau is not a ones matrix, then tau is normalized over
%      columns.
%
%OUTPUTS:
% CLL: The completed log-likelihood for the given choice of clusters. This
%      is the traditional log-likelihood if tau was empty or given as a
%      ones matrix.
% tau: A numClusters-by-numParams matrix where the i,j entry indicates
%      the jth sample originated from the ith cluster with probability
%      tau(i,j). If left empty, an all ones matrix is used which will cause
%      the returned CLL value to be the traditional log-likelihood. If
%      given but tau is not a ones matrix, then tau is normalized over
%      columns.
%
%REFERENCES:
%[1] S. Akogul and M. Erisoglu, "A Comparison of Information Criteria in
%    Clustering Based on Mixture of Multivariate Normal Distributions,"
%    Mathematical and Computational Applications, vol. 21, no. 3, p. 34,
%    Aug. 2016.
%
%August 2021 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

% Determine problem size using inputs.
numSamples = size(data,2);
numClusters = length(clusterWeights);

if numSamples==0 || numClusters==0
    CLL = NaN;
    tau = [];
    return
end

% Check that weights are normalizable.
if isrow(clusterWeights)
    clusterWeights = clusterWeights';
end
weightTotal = sum(clusterWeights);
if weightTotal~=1
    clusterWeights = clusterWeights/weightTotal;
    if any(isinf(clusterWeights))
        error('ClusterWeights must sum to 1 or be normalizable to 1.')
    end
end

% Compute the likelihoods of the data conditioned on the parameters.
density = zeros(numClusters,numSamples);
for samp = 1:numSamples
    for clust = 1:numClusters
        density(clust,samp) = likelihood(data(:,samp),params{clust,:});
    end
end

% Ensure tau is properly defined. Use Eq. 7 if not given. Use a ones matrix
% if tau is empty (traditional log-likelihood). If tau is given but is not
% a ones matrix and not left stochastic (sum over columns is 1), then tau's
% columns are normalized.
if ~exist('tau','var')
    tau = repmat(clusterWeights,1,numSamples).*density;
    tau = tau*diag(1./sum(tau,1));
elseif isempty(tau)
    tau = ones(numClusters,numSamples);
elseif ~all(size(tau)==[numClusters,numSamples])
    error('Tau must be of the shape numClusters-by-numSamples.')
elseif ~all(tau==1,'all')&&all(sum(tau,1)~=1)
    tau = tau*diag(1./sum(tau,1));
end

% Compute the conditional expected value of the completed log-likelihood
% according to Eq. 8.
tmp = tau.*(log(density)+repmat(log(clusterWeights),1,numSamples));
CLL = sum(tmp(tmp>-Inf),'all'); % Ensure invalid are not included.

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
