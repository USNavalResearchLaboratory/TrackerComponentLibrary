function [CRLB,J]=getEstimatorMinMSEBound(RInv,stateJacob,statJacob,biasVec,biasJacob,numMeasDims,checkPosDef)
%%GETESTIMATORMINMSEBOUND This evaluates a generalization of the Cramer-Rao
%       lower bound (CRLB) specifically for the case where measurements are
%       corrupted with multivariate Gaussian noise. The measurement model
%       is z=h(x)+w, where x is the (deterministic) state, h is the
%       (possible nonlinear) measurement function and w is zero-mean
%       Gaussian noise with covariance matrix R. h can vary between
%       measurements that are being fused. This offers a lower bound on the
%       mean-squared error of a biased multivariate statistic. The
%       statistic is just a transformation of the state, so zStat=f(x).
%       Biased means that if T(x) is the estimator, then the expecttation
%       E{T(x)-f(x)}=b(x), where b(x) is the bias. If one doesn't know the
%       bias and its gradient, omitting those terms provides a lower bound
%       on the covariance matrix of an unbiased statistic. If the statistic
%       Jacobian is omitted, then the bound is directly on the state x.
%       The bias term can be useful when one already has an estimator with
%       a known or approximated bias and they want to approximate its
%       accuracy with the MSE matrix bound rather than analytically
%       evaluating its true MSE.
% 
%INPUTS: RInv The zDimXzDimXnumMeas set of inverse covariance matrices
%          associated with the multivariate Gaussian noise corrupting each
%          of the numMeas measurements. If the dimensionality of the
%          measurements varies, then zDim is the maximum
%          dimensionality of any measurement and then numMeasDims is
%          required so that if sel=1:numMeasDims(k), then RInv(sel,sel,k)
%          is the submatrix used for the kth measurement.
% stateJacob The zDimXxDimXnumMeas Jacobian matrices of derivatives of the
%          measurement function h taken with respect to the elements of
%          the target state for every measurement. If a single zDimXxDim
%          matrix is passed, then it is assumed that this is the same for
%          all numMeas measurement. If the dimensionality of the
%          measurements varies, then zDim is the maximum
%          dimensionality of any measurement and then numMeasDims is
%          required so that if sel=1:numMeasDims(k), then
%          stateJacob(sel,:,k) is the submatrix used for the kth
%          measurement.
% statJacob The statDimXxDim matrix of derivatives of the the statistic
%          function f(x) taken with respect to the state x. If this is
%          omitted or an empty matrix is passed, then an xDimXxDim identity
%          matrix is used, which is the same as h(x)=x.
%  biasVec A statDimX1 vector, if provided. This is the bias of the assumed
%          estimator of the statistics. If this value and the next value
%          are omitted or empty matrices are provided, then the CRLB for an
%          unbiased statistic is computed.
% biasJacob A statDimXxDim matrix of derivatives of the estimator bias b(x)
%          taken with respect to the state x. If this is provided, then
%          biasVec must be provided. If omitted or an empty matrix is
%          passed, then this is assumed to be 0.
% numMeasDims If all measurement shave the same dimensionality , then this
%          input should be omitted or an empty matrix passed. Otherwise,
%          this is a length numMeas vector that specified how many
%          dimensions each measurement has.
% checkPosDef If this is true, then a check is performed as to whether the
%          Fisher information matrix in the algorithm is full rank.
%          If it is not, then the CRLB matrix will be returned empty and no
%          attempt will be made to invert it. The default if omitted or an
%          empty matrix is passed is false.
%
%OUTPUTS: CRLB The statDimXstatDim CRLB matrix or, if the FIM was singular,
%              and checkPosDef is true, an empty matrix.
%            J The Fisher infromation matrix that went into computing CRLB.
%
%A simple derivation of the standard multivariate CRLB is derived in [1].
%If in the derivation there, one replaces the difference T_k(x)-x_k (the
%difference of the best estimator of the kth component of x and the state)
%with T_k(x)-f_k(x), where f_k(x) is the kth component of a desired
%statistic and thus T_k(x) becomes an estimator of a component of a
%statistic, then one will derive the CRLB of a statistic instead of just of
%the state. In the scalar case, this was already in Rao's original paper in
%[2].
%
%To consider a biased statistic, one then replaces the assumption in the
%paper of a zero expected value: E{T_k(x)-x_k}=0 with
%E{T_k(x)-f_k(x)}=b_k(x). The rest of the derivation in [1] is pretty much
%the same, making sure that one is evaluating Cov{Z} and not just E{Z*Z'}.
%Completing the derivation just using E{Z*Z'} and not Cov{Z} will lead one
%to have a lower bound that is smaller by a bias*bias' term. 
%
%Note that is is possible for a biased estimator to have a lower MSE than
%an unbiased estimator. One such example is given in [3].
%
%EXAMPLE 1:
%This is a simple example, where z=x+w with w having the identity matrix as
%the covariance matrix. We take x to be 2X1. We want to estimate
%f(x)=x'*x. The assumed biased estimator is T(z)=z'*z+sum(z)/100. (Note
%that even T(z)=z'*z is biased, but we choose to make the bias non-
%constant). The bias in 2D can be found to be b(x)=2+sum(x)/100.
%The analytic mean squared error in 2D is
%MSE=(80002+40001*x(1)^2+2*x(1)*(400+x(2))+x(2)*(800+40001*x(2)))/10000
%We plot the MSE and the CRLB for a fixed x(1), while varying x(2). One can
%see a reasonably good level of agreement, though the bound is not tight.
% MSE=@(x)(80002+40001*x(1)^2+2*x(1)*(400+x(2))+x(2)*(800+40001*x(2)))/10000;
% biasVec=@(x)(2+sum(x)/100);
% biasJacob=@(x)[1/100,1/100];
% statJacob=@(x)2*x';
% stateJacob=eye(2,2);
% RInv=eye(2,2);
% numPts=200;
% x1=-1/2;
% x2=linspace(-5,5,numPts);
% MSEVal=zeros(1,numPts);
% CRLB=zeros(1,numPts);
% for curPt=1:numPts
%     xCur=[x1;x2(curPt)];
%     MSEVal(curPt)=MSE(xCur);
%     CRLB(curPt)=getEstimatorMinMSEBound(RInv,stateJacob,statJacob(xCur),biasVec(xCur),biasJacob(xCur));
% end
% figure(1)
% clf
% hold on
% plot(x2,MSEVal,'-b','linewidth',2)
% plot(x2,CRLB,'--k','linewidth',2)
% legend('Actual MSE','CRLB','location','north')
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%[2] C. R. Rao, "Information and the accuray attainable in the estimation
%    of statistical parameters," Bulletin of the Calcutta Mathematical
%    Society, vol. 37, no. 3, pp. 81-91, 1945.
%[3] P. Stoica and R. L. Moses, "On biased estimators and the unbiased
%    Cramér-Rao lower bound," Signal Processing, vol. 21, no. 4, pp. 349-
%    350, Dec. 1990.
%
%November 2022 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=size(RInv,3);
xDim=size(stateJacob,2);
statDim=size(statJacob,1);

if(nargin<7||isempty(checkPosDef))
    checkPosDef=false;
end

if(nargin<6)
    numMeasDims=[];
end

if(nargin<5||isempty(biasJacob))
    biasJacob=zeros(statDim,xDim);
end

if(nargin<4||isempty(biasVec))
    biasVec=zeros(statDim,1);
end

if(nargin<3||isempty(statJacob))
    statJacob=eye(xDim,xDim);
end

if(numMeas>1&&size(stateJacob,3)==1)
    stateJacob=repmat(stateJacob,[1,1,numMeas]);
end

%First get the FIM
J=zeros(xDim,xDim);
if(isempty(numMeasDims))
    %If all measurements are zDim in size
    for k=1:numMeas
        J=J+stateJacob(:,:,k)'*RInv(:,:,k)*stateJacob(:,:,k);
    end
else
    for k=1:numMeas
        zDimCur=numMeasDims(k);
        sel=1:zDimCur;
        J=J+stateJacob(sel,:,k)'*RInv(sel,sel,k)*stateJacob(sel,:,k);
    end
end

if(checkPosDef)
    %Check whether or not J is positive definite.
    if(matrixRank(J)<xDim)
        CRLB=[];
        return;
    end
end

CRLB=(statJacob+biasJacob)*inv(J)*(statJacob+biasJacob)'+biasVec*biasVec';

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
