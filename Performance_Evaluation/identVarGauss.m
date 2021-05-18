function [IDVarNorm,IDVar,IDVarMax]=identVarGauss(w,mu,Sigma,xDim)
%%IDENTVARGAUSS Given a Gaussian mixture representing uncertain target
%               states, determine the normalized and regular identity (ID)
%               variances, as defined in [1]. The identity variance is a
%               measure of uncertainty in the identities of which target is
%               which given the fact that the exact positions of the
%               targets are uncertain (Gaussian). The identity
%               variance provides a metric of how confident one is in
%               which target is which. This function provides an explicit
%               solution under the assumption that the joint probability
%               distribution function (PDF) of the targets is a Gaussian
%               mixture.
%
%INPUTS: w A numHypX1 or 1XnumHyp vector of weights of the hypotheses. Note
%          that w(i)>=0 for all i and sum(w)=1.
%       mu This is an (xDim*numTar)XnumHyp set of stacked state vectors for
%          all of the targets for each hypothesis. Alternatively, if the
%          target states are all uncorrelated, this can be an
%          xDimXnumTarXnumHyp hypermatrix.
%    Sigma This is an (xDim*numTar)XxDim*numTar)XnumHyp set of covariance
%          matrices for all of the stacked mean values in mu for each
%          hypothesis. If all of the target states are uncorrelated, an
%          xDimXxDimXnumTarXnumHyp set of hypermatrices for each target and
%          hypothesis individually can be passed.
%     xDim If the targets are correlated (mu is (xDim*numTar)XnumHyp in
%          size), then the state dimensions size xDim must be explicitly
%          provided.
%
%OUTPUTS: IDVarNorm The normalized ID variance. This is a value between 0
%               and 1. Zero means that the identities of the targets are
%               completely unknown, and one means that they are completely
%               certain.
%         IDVar The non-normalized ID variance.
%      IDVarMax The normalizing constant for the ID variance.
%
%This implements the algorithm of [1].
%
%EXAMPLE:
%Two targets, two hypotheses, three dimensional states.
% mu=zeros(3,2,2);
% mu(:,1,1)=[20;-30;0;];
% mu(:,1,2)=[-15;20;1];
% mu(:,2,1)=[-15;20;3];
% mu(:,2,2)=[-15;20;3];
% Sigma(:,:,1,1)=4*eye(3);
% Sigma(:,:,1,2)=Sigma(:,:,1,1);
% Sigma(:,:,2,1)=[1,  0.5, -0.5;
%                 0.5,  2,  0.5;
%                -0.5,0.5,  3];
% Sigma(:,:,2,2)=2*Sigma(:,:,2,1);
% w=[0.5;0.5];
% IDVarNorm0=identVarGauss(w,mu,Sigma)
% %One will get an ID variance of about 0.8161. However, if all first target
% %hypothese are moved far away from the second target, then the ambiguity is
% %reduced.
% mu(:,1,1)=mu(:,1,1)+500;
% mu(:,1,2)=mu(:,1,2)+500;
% IDVarNorm1=identVarGauss(w,mu,Sigma)
% %Here, one gets an ID variance of essentially 1, meaning that the
% %identities are very clear. on the other hand, if the targets are made to
% %coincide, then the identity variance becomes zero.
% mu(:,2,1)=mu(:,1,1);
% mu(:,2,2)=mu(:,1,2);
% Sigma(:,:,2,1)=Sigma(:,:,1,1);
% Sigma(:,:,2,2)=Sigma(:,:,1,2);
% IDVarNorm2=identVarGauss(w,mu,Sigma)
% %Now, the identity variance is zero.
%
%REFERENCES:
%[1] D. F. Crouse and P. Willett, "Identity variance for multi-object
%    estimation," in Proceedings of SPIE: Signal and Data Processing of
%    Small Targets, vol. 8137, San Diego, CA, 21 Aug. 2011.
%
%November 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numHyp=length(w);

if(nargin<4||isempty(xDim))
    xDim=size(mu,1);
end

if(ndims(mu)==3)
    numTar=size(mu,2);
    totalDim=xDim*numTar;
    mu=reshape(mu,totalDim,numHyp);
    
    SigmaNew=zeros(xDim*numTar,xDim*numTar,numHyp);
    for curHyp=1:numHyp
        span=1:xDim;
        for curTar=1:numTar
            SigmaNew(span,span,curHyp)=Sigma(:,:,curTar,curHyp);
            span=span+xDim;
        end
    end
    Sigma=SigmaNew;
else
    totalDim=size(mu,1);
    numTar=totalDim/xDim;
end

%The sizes of the adjusted inputs:
%w is numHypX1 or 1XnumHyp
%mu is xDimXnumTarXnumHyp
%Sigma is (xDim*numTar)X(xDim*numTar)XnumHyp

w=w(:);

%The total number of permutations.
numTarPerm=factorial(numTar);

%Here, we implement Equation 22. Note that the matrix H in Equation 31c is
%the same if i and j are swapped.
val1=0;%For the value of the first sum in Equation 22.
val2=0;%For the value of the second sum in 22 and the value of Equation 23.
for curI=1:numTarPerm
    curTerm=w'*calcHTilde(w,mu,Sigma,numTar,xDim,curI,curI)*w;
    val1=val1+curTerm;
    val2=val2+curTerm;
    
    for curJ=curI+1:numTarPerm
        curTerm=w'*calcHTilde(w,mu,Sigma,numTar,xDim,curI,curJ)*w;
        val2=val2+2*curTerm;%The 2 is for the ordering i,j as well as j,i
    end
end

IDVar=val1/numTarPerm-val2/numTarPerm^2;
IDVarMax=val2*(numTarPerm-1)/(numTarPerm^2);
IDVarNorm=IDVar/IDVarMax;
end

function Ht=calcHTilde(w,mu,Sigma,numTar,xDim,i,j)
%This function implements Equation 31c.

    numHyp=length(w);
    Ht=zeros(numHyp,numHyp);
    
    idxI=getPermIndices(i-1,numTar,xDim);
    idxJ=getPermIndices(j-1,numTar,xDim);
    
    for m=1:numHyp
        mumi=mu(idxI,m);
        Sigmami=Sigma(idxI,idxI,m);
        
        for n=m:numHyp
            munj=mu(idxJ,n);
            Sigmanj=Sigma(idxJ,idxJ,n);

            SigmamiInv=inv(Sigmami);
            SigmanjInv=inv(Sigmanj);
            
            Sigmamnij=inv(SigmamiInv+SigmanjInv);
            mumnij=Sigmamnij*(SigmamiInv*mumi+SigmanjInv*munj);

            val=0;
            for k=1:numHyp
                val=val+w(k)*GaussianD.PDF(mumnij,mu(:,k),Sigmamnij+Sigma(:,:,k));
            end
            val=val*w(m)*w(n)*GaussianD.PDF(mumi,munj,Sigmami+Sigmanj);
            
            Ht(m,n)=val;
            Ht(n,m)=val;
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
