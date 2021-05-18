function [wRed,muRed,PRed]=WestGaussReduction(w,mu,P,K,distMeas,algorithm,gammaVal,KMax)
%%WESTGAUSSREDUCTION Perform Gaussian mixture reduction as described by
%               West in Section 2.3 of [1] or via the enhanced West
%               algorithm, which is described in Section 3.2 of [2].
%
%INPUTS: w An NX1 or 1XN vector of weights of the components of the
%          original Gaussian mixture.
%       mu An xDimXN matrix of the means of the vector components of the
%          original Gaussian mixture.
%        P An xDim XxDim XN hypermatrix of the covariance matrices for the
%          components of the original Gaussian mixture.
%        K The number of components desired in the mixture after reduction.
% distMeas The algorithm sequentially merges the Gaussian mixture component
%          having the least weight (or modified weight in the enhanced West
%          algorithm) with the nearest other component. This specifies the
%          distance measure to use to determine what is "near." Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            Kullback-Leiber (KL) divergence.
%          1 Use the integrated squared error (ISE).
% algorithm An optional parametering indicating the algorithm to use.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            West algorithm of Section 2.3 of [1].
%          1 Use the enhanced west of [2]. This the same same as the West
%            algorithm, except the weights used for choosing which
%            component to merge are a little different.
% gammaVal This is an optional threshold introduced in [2]. This is maximum
%          cost allowed for two components to be merged. If at any time,
%          the minimum cost is larger than this, the mixture redution stops
%          early with possibly more than K components. If this parameter is
%          omitted or an empty matrix is passed, then reduction is always
%          done until K components is reached.
%     KMax When gammaVal is given, there can be a possibility that too many
%          components are retained. This is the maximum number of
%          components to allow. Stopping due to gammaVal will not occur
%          unless the number of components left is less than or equal to
%          KMax. The default if omitted or an empty matrix is passed is Inf
%          (no limit).
%
%OUTPUTS: wRed The KX1 weights of the mixture after reduction. In this and
%              the other outputs, there can me more than K components if
%              the 
%        muRed The xDimXK means of the mixture after reduction.
%         PRed The xDimXxDimXK covariance matrices of the mixture after
%              reduction.
%
%EXAMPLE:
%This scalar example plots a full 9-component PDF and then the PDF reduced
%to 6 components using sequential brute-force reudction, West's algorithm
%and also Runnals' algorithm.
% w=[0.03,0.18,0.12,0.19,0.16,0.06,0.1,0.08,0.06];
% w=w/sum(w);
% n=length(w);
% mu=[1.45,2.20,0.67,0.48,0.91,1.01,1.42,2.77,0.89];
% P=[0.0487,0.0305,0.1171,0.0174,0.0102, 0.0323, 0.0380, 0.0115, 0.0679];
% P=reshape(P,[1,1,n]);
% 
% k=6;%Number of reduced components.
% %Sequential brute-force reduction.
% [wRedBF,muRedBF,PRedBF]=bruteForceGaussMixRed(w,mu,P,k,true);
% 
% distMeas=0;
% algorithm=0;
% %West's reduction
% [wRedW,muRedW,PRedW]=WestGaussReduction(w,mu,P,k,distMeas,algorithm);
% %Runnal's reduction
% [wRedR,muRedR,PRedR]=RunnalsGaussMixRed(w,mu,P,k);
% 
% numPoints=500;
% xVals=linspace(0,3,numPoints);
% PDFVals1=GaussianMixtureD.PDF(xVals,w,mu,P);
% PDFVals2=GaussianMixtureD.PDF(xVals,wRedBF,muRedBF,PRedBF);
% PDFVals3=GaussianMixtureD.PDF(xVals,wRedW,muRedW,PRedW);
% PDFVals4=GaussianMixtureD.PDF(xVals,wRedR,muRedR,PRedR);
% figure(1)
% clf
% hold on
% plot(xVals,PDFVals1,'-k','linewidth',4)
% plot(xVals,PDFVals2,'--r','linewidth',2)
% plot(xVals,PDFVals3,'-.g','linewidth',2)
% plot(xVals,PDFVals4,'-m','linewidth',1)
% legend('Full Mixture','Sequential Brute-Force','Algorithm of West','Algorithm of Runnals')
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] M. West, "Approximate posterior distributions by mixture," Journal of
%    the Royal Statistical Society. Series B (Methodological), vol. 55, no.
%    2, pp. 409-422, 1993.
%[2] H. D. Chen, K. C. Chang, and C. Smith, "Constrained optimized weight
%    adaptation for Gaussian mixture reduction," in Proceedings of SPIE:
%    Signal Processing, Sensor Fusion, and Target Recognition XIX, vol.
%    7697, Orlando, FL, 27 Apr. 2010.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<8||isempty(KMax))
    KMax=Inf; 
end

if(nargin<7||isempty(gammaVal))
    gammaVal=Inf; 
end

if(nargin<6||isempty(algorithm))
    algorithm=0; 
end

if(nargin<5||isempty(distMeas))
    distMeas=0;
end

if(algorithm~=0&&algorithm~=1)
    error('Invalid algorithm specified.')
end

if(distMeas~=0&&distMeas~=1)
    error('Invalid distance measure specified.')
end

N=length(w);
xDim=size(mu,1);

if(N<=K)%If the mixture is already reduced.
    wRed=w;
    muRed=mu;
    PRed=P;
    return
end

%Both cost measures will require determinants.
detVec=zeros(N,1);
for i=1:N
    detVec(i)=det(P(:,:,i));
end

if(distMeas==1)
%Precompute values of Jrr in the integrated squared error (ISE) cost
%that is used so that they do not have to be recomputed every time an
%element is deleted.
    Jrr=((4*pi)^xDim*detVec).^(-1/2);
end

if(algorithm~=0)
    costVec=w./detVec;
else
    costVec=w;
end

idx=1:N;
for r=N:-1:(K+1)
    [~,minIdx1]=min(costVec(idx(1:r)));

    minIdxFull=idx(minIdx1);

    %Now, we find the nearest neighbor of the minimum weight Gaussian and 
    %merge the minimum weight Gaussian with its nearest neighbor.
    w1=w(minIdxFull);
    mu1=mu(:,minIdxFull);
    P1=P(:,:,minIdxFull);
    minDist=Inf;
    minIdx2=[];
    for i=1:r
        if(i==minIdx1)
            continue;
        end
        idxCur=idx(i);
        
        muCur=mu(:,idxCur);
        PCur=P(:,:,idxCur);
        
        if(distMeas==1)%If using the ISE.
            JrrCur=Jrr(idxCur);
            curDist=ISESimp(mu1,P1,muCur,PCur,JrrCur);
        else%If using the KL Divergence.
            detP1=detVec(minIdxFull);
            detP2=detVec(idxCur);
            
            curDist=KLDistSimp(mu1,P1,muCur,PCur,detP1,detP2);
        end
        
        if(curDist<minDist)
            minDist=curDist;
            minIdx2=i;
        end
    end
    
    if(minDist>=gammaVal&&r<=KMax)
        %Stop the reduction process; the minimum distance is larger than
        %the threshold.
        K=r;
        break;
    end

    %Merge the first Gaussian with its nearest neighbor and put the
    %result in the first index location. The second index is deleted.
    minIdxFull2=idx(minIdx2);
    w2=w(minIdxFull2);
    mu2=mu(:,minIdxFull2);
    P2=P(:,:,minIdxFull2);
    
    wMerged=w1+w2;
    lambda1=w1/wMerged;
    lambda2=w2/wMerged;
    muMerged=lambda1*mu1+lambda2*mu2;
    diffVal=mu1-mu2;
    PMerged=lambda1*P1+lambda2*P2+lambda1*lambda2*(diffVal*diffVal');

    %Overwrire the minimum weight component.
    w(minIdxFull)=wMerged;
    mu(:,minIdxFull)=muMerged;
    P(:,:,minIdxFull)=PMerged;
    
    detCur=det(PMerged);
    if(distMeas==1)
        JrrCur(minIdxFull)=((4*pi)^xDim*detCur).^(-1/2);
    else
        detVec(minIdxFull)=detCur;
    end

    if(algorithm~=0)
        costVec(minIdxFull)=wMerged/detCur;
    else
        costVec(minIdxFull)=costVec(minIdxFull)+costVec(minIdxFull2);
    end

    %Delete the nearest neighbor component.
    idx(minIdx2)=idx(r);
end

sel=idx(1:K);
wRed=w(sel);
muRed=mu(:,sel);
PRed=P(:,:,sel);

end

function ISEVal=ISESimp(mu1,P1,mu2,P2,Jrr)
%%ISESIMP The integrated squared error (ISE) between two multivariate
%         Gaussian PDFs, but if we assume that this function is called
%         multiple times for a constant mu1, P1 and varying mu2, P2, we do
%         not need the constant term involving only mu1 and P1, so it is
%         omitted here.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    Jhr=GaussianD.PDF(mu1,mu2,P1+P2);

    %The Jrr value should be equal to GaussianD.PDF(mu2,mu2,2*P2)
    ISEVal=Jrr-2*Jhr;
end

function KLDist=KLDistSimp(mu1,P1,mu2,P2,detP1,detP2)
%%KLDISTSIMP The Kullback-Leiber (KL) divergence between two multivariate
%            Gaussian PDFs. We assume that the determinants of P1 and P2
%            have been precomputed and are provided to this function.
%
%May 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.

    KLDist=0.5*(trace(lsqminnorm(P2,(P1-P2+(mu1-mu2)*(mu1-mu2)')))+log(detP2)-log(detP1));
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
