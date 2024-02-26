function [CLL,tau] = completedLogLikelihood(data,likelihood,clusterWeights,params,tau)
%%COMPLETEDLOGLIKELIHOOD Computes the conditional expected value of the
%                        completed log-likelihood (CLL) of a cluster model.
%                        This is part of the E step of the expectation
%                        maximization algorithm described in [1]. To
%                        compute the log-likelihood, give tau as an empty
%                        matrix.
%
%INPUTS: data A numDim-by-numSamples matrix where numDim is the state
%            dimension and numSamples is the number of vectors
%            representing individual data points.
% likelihood A function handle which takes as inputs a single data vector
%            from the input data and a single row of parameters from the
%            input params in that order.
% clusterWeights A vector of normalizable weights with length equal to the
%            number of clusters in the model (numClusters).
%     params A numClusters-by-numParams cell array where each row
%            corresponds to the model parameters for the corresponding
%            cluster and params{clust,:} should unpack the model parameters
%            in the correct order to be given to likelihood.
%        tau A numClusters-by-numParams matrix where the i,j entry
%            indicates the jth sample originated from the ith cluster with
%            probability tau(i,j). If left empty, an all ones matrix is
%            used which will cause the returned CLL value to be the
%            traditional log-likelihood. If given but tau is not a ones
%            matrix, then tau is normalized over columns.
%
%OUTPUTS: CLL The completed log-likelihood for the given choice of
%             clusters. This is the traditional log-likelihood if tau was
%             empty or given as a ones matrix.
%        tau  A numClusters-by-numParams matrix where the i,j entry
%             indicates the jth sample originated from the ith cluster with
%             probability tau(i,j). If left empty, an all ones matrix is
%             used which will cause the returned CLL value to be the
%             traditional log-likelihood. If given but tau is not a ones
%             matrix, then tau is normalized over columns.
%
%EXAMPLE 1: Computes the CLL for a collection of 10 samples.
% w=[0.03, 0.18, 0.12, 0.19, 0.02, 0.16, 0.06, 0.1, 0.08, 0.06];
% mu=[1.45, 2.20, 0.67, 0.48, 1.49, 0.91, 1.01, 1.42, 2.77, 0.89];
% P=reshape([0.0487, 0.0305, 0.1171, 0.0174, 0.0295,0.0102, 0.0323, 0.0380, 0.0115, 0.0679],[1,1,10]);
% [~,~,~,~,~,~,wAll,muAll,PAll]=RunnalsGaussMixRed(w,mu,P,1,Inf);
% CLL = zeros(1,length(w));
% for k = 1:length(w)
%     idx=length(w)-k+1;
%     sel=1:k;
%     wRed=wAll(sel,idx);
%     muRed=muAll(:,sel,idx);
%     PRed=PAll(:,:,sel,idx);
% 
%     params = cell(size(muRed,2),2);
%     for clust = 1:size(muRed,2)
%         params{clust,1} = muRed(:,clust);
%         params{clust,2} = PRed(:,:,clust);
%     end
%     [CLL(k),~] = completedLogLikelihood(mu,@(z,mu,R)GaussianD.PDF(z,mu,R),wRed,params);
% end
% figure(1);clf;hold on;
% plot(CLL,'b')
% xlabel('# of Merged Clusters')
% ylabel('CLL')
%
%EXAMPLE 2: Generates samples from a mixture of 3 Gaussians and computes
%           the CLL. The reduction which maximizes the CLL is then plotted
%           against the true mixture.
% mu1 = 2; P1 = 0.1;
% mu2 = 5; P2 = 0.1;
% mu3 = 10; P3 = 1;
% p = 1;
% num1 = 30;
% num2 = 30;
% num3 = 30;
% num = num1+num2+num3;
% mu = zeros(1,num);
% mu(1:num1) = GaussianD.rand(num1,mu1,P1);
% P(1,1,1:num1) = p*ones(1,1,num1);
% mu(num1+1:num1+num2) = GaussianD.rand(num2,mu2,P2);
% P(1,1,(num1+1):(num1+num2)) = p*ones(1,1,num2);
% mu(num1+num2+1:num1+num2+num3) = GaussianD.rand(num3,mu3,P3);
% P(1,1,(num1+num2+1):(num1+num2+num3)) = p*ones(1,1,num3);
% w = ones(1,num);
% [~,~,~,~,~,~,wAll,muAll,PAll]=RunnalsGaussMixRed(w,mu,P,1,Inf);
% CLL = zeros(1,length(w));
% for k = 1:length(w)
%     idx=length(w)-k+1;
%     sel=1:k;
%     wRed=wAll(sel,idx);
%     muRed=muAll(:,sel,idx);
%     PRed=PAll(:,:,sel,idx);
% 
%     params = cell(size(muRed,2),2);
%     for clust = 1:size(muRed,2)
%         params{clust,1} = muRed(:,clust);
%         params{clust,2} = PRed(:,:,clust);
%     end
%     [CLL(k),~] = completedLogLikelihood(mu,@(z,mu,R)GaussianD.PDF(z,mu,R),wRed,params);
% end
% figure(1);clf;hold on;
% plot(CLL,'b')
% xlabel('# of Merged Clusters')
% ylabel('CLL')
% figure(2);clf;hold on;
% plot(linspace(-5,15),GaussianMixtureD.PDF(linspace(-5,15),w,mu,P))
% plot(linspace(-5,15),GaussianMixtureD.PDF(linspace(-5,15),wAll(1:3,88),muAll(1,1:3,88),PAll(:,:,1:3,88)))
% legend('Truth','Reduced Samples')
% title('PDFs of Mixtures')
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
