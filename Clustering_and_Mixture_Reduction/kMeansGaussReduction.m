function [wRed,xRed,PRed,exitCode]=kMeansGaussReduction(wOrig,xOrig,POrig,param4,xInit,PInit,maxIter)
%%KMEANSGAUSSREDUCTION Perform Gaussian mixture reduction by iterating a k-
%              means algorithm given an inital reduced estimator (or using
%              Runnals' algorithm for the initial estimate if none is
%              given. 
%
%INPUTS: wOrig An NX1 or 1XN vector of weights of the components of the
%          original Gaussian mixture.
%   muOrig An xDimXN matrix of the means of the vector components of the
%          original Gaussian mixture.
%    POrig An xDimXxDim XN hypermatrix of the covariance matrices for the
%          components of the original Gaussian mixture.
%   param4 If not initial estimate of a reduced distribution is provided,
%          then this is k, the number of components to which the
%          distribution should be reduced (and xInit and PInit are either
%          not provided or are empty matrices) and RunnalsGaussMixRed will
%          be used to obtain an initial estimate. Otherwise, if an initial
%          reduced distribution is provided, this is wInit, the kX1 or 1Xk
%          set of weights for the initial estimate of the reduced
%          distribution.
%    xInit If an initial reduced distribution is provided, this is the
%          xDimXk set of Gaussian mean vectors. Otherwise, this parameter
%          can be omitted or an empty matrix passed.
%    PInit If an initial reduced distribution is provided, this is the
%          xDimXxDimXk set of covariance matrices. Otherwise, this
%          parameter can be omitted or an empty matrix passed.
%  maxIter The maximum number of iterations of the k means algorithm to
%          perform. The default if omitted or an empty matrix is passed is
%          50.
%
%OUTPUTS: wRed The KX1 weights of the mixture after reduction.
%        muRed The xDimXK means of the mixture after reduction.
%         PRed The xDimXxDimXK covariance matrices of the mixture after
%              reduction.
%     exitCode A parameter indicating how the algorithm terminated.
%              Possible values are:
%              0 The algorithm converged.
%              1 the maximum number of iterations elapsed.
%
%This implements the algorithm of [1] and [2], which is also described in
%[3]. TGlobal convergence is not guaranteed and the performance depends a
%lot on the quality of the initial estiamte. The Kullback-Leiblber
%divergence is used as a closeness criterion between components via the
%KLDivGauss function. Additional algorithmic optimizations precomputing
%various determinants are possible beyond what was done here.
%
%EXAMPLE:
%We reduce a 10-component mixture to five components. We first plot
%Runnals' algorithm. We then use a perturbed version of the Runnals
%solution as the initial estimate to the k-means estimator and we find that
%it converges back to the Runnals solution.
% w=[0.03,0.18,0.12,0.19,0.02,0.16,0.06,0.1,0.08,0.06];
% n=length(w);
% mu=[1.45,2.20,0.67,0.48,1.49,0.91,1.01,1.42,2.77,0.89];
% P=[0.0487,0.0305,0.1171,0.0174,0.0295,0.0102, 0.0323, 0.0380, 0.0115, 0.0679];
% P=reshape(P,[1,1,n]);
% 
% k=6;%Number of reduced components.
% 
% %Runnal's reduction
% [wRedR,muRedR,PRedR]=RunnalsGaussMixRed(w,mu,P,k);
% 
% %Perturb Runnal's solution so that the kMeansReduction has a bad
% %initialization.
% wRedRP=wRedR;
% wRedRP(1)=wRedR(1)+0.1;
% wRedRP=wRedRP/sum(wRedRP);
% muRedRP=muRedR+0.1;
% [wRedkM,muRedkM,PRedkM]=kMeansGaussReduction(w,mu,P,wRedRP,muRedRP,PRedR);
% 
% numPoints=500;
% xVals=linspace(0,3,numPoints);
% PDFVals1=GaussianMixtureD.PDF(xVals,w,mu,P);
% PDFVals2=GaussianMixtureD.PDF(xVals,wRedR,muRedR,PRedR);
% PDFVals3=GaussianMixtureD.PDF(xVals,wRedRP,muRedRP,PRedR);
% PDFVals4=GaussianMixtureD.PDF(xVals,wRedkM,muRedkM,PRedkM);
% figure(1)
% clf
% hold on
% plot(xVals,PDFVals1,'-k','linewidth',4)
% plot(xVals,PDFVals2,'-r','linewidth',4)
% plot(xVals,PDFVals3,'-.c','linewidth',2)
% plot(xVals,PDFVals4,'--g','linewidth',2)
% legend('Full Mixture','Algorithm of Runnals','Initialization of k-Means','k-Means Algorithm')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] D. Schieferdecker and M. F. Huber, "Gaussian mixture reduction via
%    clustering," in Proceedings of the 12th International Conference on
%    Information Fusion, Seattle, WA, 6-9 Jul. 2009, pp. 1536-1543.
%[2] A. Nikseresht and M. Gelgon, "Gossip-based computation of a Gaussian
%    mixture model for distributed multimedia indexing," IEEE Transactions
%    on Multimedia, vol. 10, no. 3, pp. 385-392, Apr. 2008.
%[3] D. F. Crouse, P. Willett, K. Pattipati, and L. Svensson, "A look at
%    Gaussian mixture reduction algorithms," in Proceedings of the 14th
%    International Conference on Information Fusion, Chicago, IL, 5-8 Jul.
%    2011.
%
%May 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(maxIter))
   maxIter=20; 
end

if(nargin<5||isempty(xInit))
%If no initial estimate is given but the desired number of components is
%given, then get an initial estimate using Runnall's algorithm.
    [wInit,xInit,PInit]=RunnalsGaussMixRed(wOrig,xOrig,POrig,param4);
else
    wInit=param4;
end

%The original number of mixture parameters.
N=length(wOrig);
%The reduced number of mixture parameters.
k=length(wInit);

cAssign=zeros(N,1);
cAssignNew=zeros(N,1);

%We assign the Gaussians in the full mixture to cluster centers in the
%current reduced mixture. Then, we recompute the cluster centers. These two
%steps are looped until the assignments no longer change or a maximum
%number of iterations has passed.

wRed=wInit;
xRed=xInit;
PRed=PInit;

exitCode=1;
for curIter=1:maxIter
    %For each point, find the closest center.
    for curP=1:N
        minCost=Inf;
        for curC=1:k
            cost=KLDistGauss(xOrig(:,curP),POrig(:,:,curP),xRed(:,curC),PRed(:,:,curC));
            if(cost<minCost)
                minCost=cost;
                cAssignNew(curP)=curC;
            end
        end
    end

    %Terminate when the assignments no longer change.
    if(sum(cAssignNew~=cAssign)==0)
        exitCode=0;%Convergence attained.
        break;
    end
    cAssign=cAssignNew;

    %Now, merge all points having common centers.
    [cSorted,idx]=sort(cAssign,'ascend');
    [~,numReps]=runLenEncode(cSorted);
    startIdx=1;
    for curC=1:k
        selIdx=idx(startIdx:(startIdx+numReps(curC)-1));
        
        wSelNorm=wOrig(selIdx);
        wSum=sum(wSelNorm);
        wSelNorm=wSelNorm/wSum;
        
        [xMerged,PMerged]=calcMixtureMoments(xOrig(:,selIdx),wSelNorm,POrig(:,:,selIdx));

        wRed(curC)=wSum;
        xRed(:,curC)=xMerged;
        PRed(:,:,curC)=PMerged;
        
        startIdx=startIdx+numReps(curC);
    end

    wRed=wRed/sum(wRed);
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
