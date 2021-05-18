classdef GaussianMixtureD
%%GAUSSIANMIXTURED Function to handle scalar or multivariate real Gaussian
%                  mixture distributions.
%Implemented methods are: mean, cov, PDF, PDFS, PDFI, rand, randS.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
methods(Static)
function meanVal=mean(w,mu)
%%MEAN Obtain the mean of the Gaussiax mixture distribution.
%
%INPUTS: w The 1XN or NX1 set of weights of the distribution components
%          such that all w>=0 and sum(w)=1.
%       mu The xDimXN set of mean values of the distribution
%          components.
%
%OUTPUTS: meanVal The xDimX1 mean of the Gaussian mixture distribution.
%
%It is simple to show that the mean of the distribution is the weighted
%sum of the means of the components.
%
%EXAMPLE:
%In this example, 
% w=[0.4,0.6];
% mu=[1,-0.5;
%    -1, 1];
% Sigma=zeros(2,2,2);
% Sigma(:,:,1)=[4/9,  14/45;
%               14/45,4/9];
% Sigma(:,:,2)=[4/9, 0;
%                 0, 4/9];
% N=10000;
% analyticMean=GaussianMixtureD.mean(w,mu);
% sampMean=mean(GaussianMixtureD.rand(N,w,mu,Sigma),2);
% RelErr=abs((sampMean-analyticMean)./analyticMean)
%In this example, the relative error will typically be 3% or better.
%Increasing the number of samples will typically decrease the relative
%error, indicating that the analytic and sample means are converging.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    meanVal=sum(bsxfun(@times,w(:).',mu),2);
end

function [covMat,meanVal]=cov(w,mu,Sigma)
%%COV Obtain the covariance matrix of the Gaussian mixture distribution
%     (the variance if it is scalar).
%
%INPUTS: w The 1XN or NX1 set of weights of the distribution components
%          such that all w>=0 and sum(w)=1.
%       mu The xDimXN set of mean values of the distribution
%          components.
%    Sigma The xDimXxDimXN set of covariance matrices of the
%          distribution components.
%
%OUTPUTS: covMat The xDimXxDim covariance matrix of the Gaussian
%                mixture.
%        meanVal The xDimX1 mean value of the distribution.
%
%The function calcMixtureMoments is called to find the covariance
%matrix.
%
%EXAMPLE:
% w=[0.4,0.6];
% mu=[1,-0.5;
%    -1, 1];
% Sigma=zeros(2,2,2);
% Sigma(:,:,1)=[4/9,  14/45;
%               14/45,4/9];
% Sigma(:,:,2)=[4/9, 0;
%                 0, 4/9];
% N=10000;
% analyticCov=GaussianMixtureD.cov(w,mu,Sigma);
% [~,sampleCov]=calcMixtureMoments(GaussianMixtureD.rand(N,w,mu,Sigma));
% RelErr=abs((sampleCov-analyticCov)./analyticCov)
%The relative error of each element will tend to be less than 1% and
%will decrease as the number of samples increases
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    [meanVal,covMat]=calcMixtureMoments(mu,w,Sigma);
end

function vals=PDF(z,w,mu,Sigma)
%%PDF Evaluate the PDF of a scalar or multivariate Gaussian
%     distribution at specified points.
%
%INPUTS: z The xDimXnumPoints set of points at which the PDF should be
%          evaluated.
%        w The 1XN or NX1 set of weights of the distribution components
%          such that all w>=0 and sum(w)=1.
%       mu The xDimXN set of mean values of the distribution
%          components.
%    Sigma The xDimXxDimXN set of covariance matrices of the
%          distribution components.
%
%OUTPUTS: vals The 1XN set of PDF values of the Gaussian mixture PDF
%              evaluated at the points in z.
%
%The PDF is the weighted sum of the PDFs of the individual components.
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing
%the PDF plot with a histogram of the random samples.
% w=[0.4,0.6];
% mu=[1,-0.5];
% sigma=zeros(1,1,2);
% sigma(:,:,1)=4/9;
% sigma(:,:,2)=5/9;
% N=5000;
% figure(1)
% clf
% histogram(GaussianMixtureD.randS(N,w,mu,sigma),'Normalization','pdf','BinLimits',[-3,3])
% hold on
% numPoints=1000;
% x=linspace(-3,3,numPoints);
% vals=GaussianMixtureD.PDF(x,w,mu,sigma.^2);
% plot(x,vals,'linewidth',2)
% axis([-3,3,0,0.45])
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    n=length(w);
    numMeas=size(z,2);

    vals=zeros(1,numMeas);
    for k=1:n
        vals=vals+w(k)*GaussianD.PDF(z,mu(:,k),Sigma(:,:,k));
    end
end

function vals=PDFS(z,w,mu,S)
%%PDF Evaluate the PDF of a scalar or multivariate Gaussian
%     distribution at specified points given the weights, means and
%     lower-triangular covariance matrices of the individual
%     components.
%
%INPUTS: z The xDimXnumPoints set of points at which the PDF should be
%          evaluated.
%        w The 1XN or NX1 set of weights of the distribution components
%          such that all w>=0 and sum(w)=1.
%       mu The xDimXN set of mean values of the distribution
%          components.
%        S The xDimXxDimXN set of lower-triangular square-root
%          covariance matrices of the distribution components.
%
%OUTPUTS: vals The 1XN set of PDF values of the Gaussian mixture PDF
%              evaluated at the points in z.
%
%The PDF is the weighted sum of the PDFs of the individual components.
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing
%the PDF plot with a histogram of the random samples.
% w=[0.4,0.6];
% mu=[1,-0.5];
% sigma=zeros(1,1,2);
% sigma(:,:,1)=4/9;
% sigma(:,:,2)=5/9;
% N=5000;
% figure(1)
% clf
% histogram(GaussianMixtureD.randS(N,w,mu,sigma),'Normalization','pdf','BinLimits',[-3,3])
% hold on
% numPoints=1000;
% x=linspace(-3,3,numPoints);
% vals=GaussianMixtureD.PDFS(x,w,mu,sigma);
% plot(x,vals,'linewidth',2)
% axis([-3,3,0,0.45])
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    n=length(w);
    numMeas=size(z,2);

    vals=zeros(1,numMeas);
    for k=1:n
        vals=vals+w(k)*GaussianD.PDFS(z,mu(:,k),S(:,:,k));
    end
end

function vals=PDFI(z,w,mu,SigmaInv,SigmaInvDet)
%%PDFI Evaluate the PDF of a scalar or multivariate Gaussian
%      distribution at specified points given the weights, means and
%      the inverses of the covariance matrices of the individual
%      components.
%
%INPUTS: z The xDimXnumPoints set of points at which the PDF should be
%          evaluated.
%        w The 1XN or NX1 set of weights of the distribution components
%          such that all w>=0 and sum(w)=1.
%       mu The xDimXN set of mean values of the distribution
%          components.
% SigmaInv The xDimXxDimXN set of inverse of covariance matrices of the
%          distribution components.
% SigmaInvDet Optionally, a length-N set of determinants of the matrices
%          in SigmaInv can be passed so as to speed up the computation. If
%          omitted or an empty matrix is passed, determinants will be taken
%          as needed.
%
%OUTPUTS: vals The 1XN set of PDF values of the Gaussian mixture PDF
%              evaluated at the points in z.
%
%The PDF is the weighted sum of the PDFs of the individual components.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    n=length(w);
    numMeas=size(z,2);

    vals=zeros(1,numMeas);
    if(nargin<5||isempty(SigmaInvDet))
        for k=1:n
            vals=vals+w(k)*GaussianD.PDFI(z,mu(:,k),SigmaInv(:,:,k));
        end
    else
        for k=1:n
            vals=vals+w(k)*GaussianD.PDFI(z,mu(:,k),SigmaInv(:,:,k),SigmaInvDet(k));
        end
    end
end

function vals=CDF(z,w,mu,varVals)
%%CDF Evaluate the cumulative distribution function (CDF) of a scalar
%     Gaussian mixture at the specified points.
%
%INPUTS: z A matrix of the point(s) at which the CDF should be
%          evaluated.
%        w The 1XN or NX1 set of weights of the distribution components
%          such that all w>=0 and sum(w)=1.
%       mu The length N set of mean values of the distribution
%          components.
%  varVals The length-N set of variances of the distribution
%          components.
%
%OUTPUTS: val The scalar value(s) of the Gaussian mixture distribution
%             evaluated at the point(s) z.
%
%The CDF of the mixture is the weighted sum of the CDFs of the
%individual Gaussian components.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    n=length(w);

    vals=zeros(size(z));
    for k=1:n
        vals=vals+w(k)*GaussianD.CDF(z,mu(k),varVals(k));
    end
end

function vals=rand(N,w,mu,P)
%%RAND Generate multivariate Gaussian mixture samples with given
%      weights, means and lower-triangular square root covariance
%      matrices.
%    
%INPUTS: N The number of random variables to generate.
%        w The 1XN or NX1 set of weights of the distribution components
%          such that all w>=0 and sum(w)=1.
%       mu The xDimXN set of mean values of the distribution
%          components.
%        P The xDimXxDimXN set of covariance matrices of the
%          distribution components. If all of them are the same, then a
%          single xDimXxDim matrix can be passed.
%
%OUTPUT: x An xDimXN matrix of random instances of the multivariate
%          Gaussian distribution.
%
%The distribution is sampled by sampling the empirical distribution of
%weights, to determine which Gaussian component to sample, and then
%sampling the chosen Gausian distribution. Lower-triangular square
%roots of the componets must be found for sampling, so if these are
%available, then the randS method will be faster.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %The number of components of the mixture.
    n=length(w);
    xDim=size(mu,1);

    if(size(P,3)==1)
        P=repmat(P,[1,1,n]); 
    end
    
    vals=zeros(xDim,N);

    %Determine which components are to be sampled.
    idx=EmpiricalD.rand([N,1],1:n,w);

    %Sample the Gaussian components that were selected.
    S=zeros(xDim,xDim,n);
    componentUsed=false(n,1);
    for curSamp=1:N
        curComp=idx(curSamp);
        %Get the lower-triangular square roots only for the components
        %that are used.
        if(componentUsed(curComp)==false)
            S(:,:,curComp)=chol(P(:,:,curComp),'lower');
            componentUsed(curComp)=true;
        end

        vals(:,curSamp)=mu(:,curComp)+S(:,:,curComp)*randn(xDim,1);
    end
end

function vals=randS(N,w,mu,S)
%%RANDS Generate multivariate Gaussian mixture samples with given
%       weights, means and lower-triangular square root covariance
%       matrices.
%
%INPUTS: N The number of random variables to generate.
%        w The 1XN or NX1 set of weights of the distribution components
%          such that all w>=0 and sum(w)=1.
%       mu The xDimXN set of mean values of the distribution
%          components.
%        S The xDimXxDimXN set of lower-triangular square-root
%          covariance matrices of the distribution components.
%
%OUTPUT: x An xDimXN matrix of random instances of the multivariate
%          Gaussian distribution.
%
%The distribution is sampled by sampling the empirical distribution of
%weights, to determine which Gaussian component to sample, and then
%sampling the chosen Gausian distribution.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %The number of components of the mixture.
    n=length(w);
    xDim=size(mu,1);

    if(size(S,3)==1)
        S=repmat(S,[1,1,n]); 
    end
    
    vals=zeros(xDim,N);

    %Determine which components are to be sampled.
    idx=EmpiricalD.rand([N,1],1:n,w);

    %Sample the Gaussian components that were selected.
    for curSamp=1:N
        curComp=idx(curSamp);

        vals(:,curSamp)=mu(:,curComp)+S(:,:,curComp)*randn(xDim,1);
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
