function [xEst,PML]=estLinShrinkageMean(z,H,R,scaleUnknown)
%%ESTLINSHRINKAGEMEAN Estimate the unknown mean of a batch of linear
%                  measurements of a multivariate state with additive
%                  Gaussian noise using an estimator that beats the
%                  expected value in mean-squared error if the
%                  dimensionality of the state is greater than 2. For a
%                  state that is two or fewer dimensions, the estimate is
%                  the expected value (the ML estimate). For more
%                  dimensions, a James-Stein (shrinkage) estimator is used.
%                  The James-Stein estimator has a lower MSE then the ML
%                  estimate, but it is biased. The measurement model is
%                  z=H*x+w where z is all of the measurements, stacked, H
%                  is all of the measurement matrices, stacked, and x is
%                  the unknown state that is to be estimated. w is additive
%                  Gaussian noise with covariance matrix R. R only needs to
%                  be known within a scalar constant of the true value.
%
%INPUTS: z A zDim X 1 vector of a measurement or a stacked vector for
%          multiple measurements.
%        H A zDim X xDim measurement matrix where zDim>=xDim. H'*R*H' 
%          should be invertible.
%        R The covariance matrix for z. This only needs to be known within
%          a scalar constant of the true value.
%scaleUnknown A boolean value. This is true if the covariance matrix R is
%          only known within a constant of the true value, otherwise, this
%          is false. To use scaleUnknown=true, it is required that
%          zDim>xDim. The default if this parameter is omitted is false.
%
%OUTPUTS: xEst The James-Stein estimator of the unknown mean. If xDim<=2,
%             this will be the ML estimate. Otherwise, this will have a MSE
%             less than the ML estimate but be biased.
%         PML The covariance matrix of the ML estimate. This can be used as
%             an upper bound on the MSE of xJS. If R is only known within a
%             scale factor, then PML is only known within the same scale
%             factor.
%
%If the scaling of R is unknown, that is if one only has R/sigma^2, where
%sigma is some unknown scale factor, then one should use scaleUnknown=true.
%When the dimensionality of x is 2 or less, the unknown scale factor does
%not matter.
%
%The algorithm is implemented from the beginning of [1] and the equation
%numbers in the comments refer to equations in that paper.
%
%Many do not realize that estimators with lower mean-squared error than the
%expected value can exist. The first proof of this fact was given in [2]
%and surprised many statisticians. The name "James-Stein" for the
%estimation algorithm is due to the paper [3].
%
%Note that when just using random measurement and measurement covariance
%matrices H and R, the algorithm will be the ML estimate, since the
%criteria for selecting a shrinkage factor require that. Code to simulate a
%scenario where the James-Stein estimator is different and outperforms the
%ML estimator is
% xDim=60;
% zDim=3*xDim;
% numRuns=5000;
% H=repmat(eye(xDim),[3,1]);
% R=eye(zDim);
% 
% SErrJS=0;
% SErrML=0;
% for curRun=1:numRuns
%     z=H*x+chol(R,'lower')*randn(zDim,1);
%     %The ML estimator
%     xML=(H'*inv(R)*H)\H'*(R\z);
%     diff=x-xML;
%     SErrML=SErrML+diff'*diff;
%     
%     xJS=estLinShrinkageMean(z,H,R,false);
%     diff=x-xJS;
%     SErrJS=SErrJS+diff'*diff;
% end
% %Display the mean squared errors under each algorithm.
% SErrJS=SErrJS/numRuns
% SErrML=SErrML/numRuns
%
%Typical simulated squared errors are
%SErrJS=20.072491066769956
%and
%SErrML=20.077663028980357
%As can be seen, the difference is often small.
%
%REFERENCES:
%[1] J. H. Manton, V. Krishnamurthy, and H. V. Poor, "James-Stein state
%    filtering algorithms," IEEE Transactions on Signal Processing, vol.
%    46, no. 9, pp. 2431-2447, Sep. 1998.
%[2] C. Stein, "Inadmissibility of the usual estimator for the mean of a
%    multivariate normal distribution," in Proceedings of the Third
%    Berkeley Symposium on Mathematical Statistics and Probability, vol. I,
%    1956, pp. 197-206.
%[3] W. James and C. Stein, "Estimation with quadratic loss," in
%    Proceedings of the Fourth Berkeley Symposium on Mathematical
%    Statistics and Probability, vol. I, 1961, pp. 361-379.
%
%May 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    scaleUnknown=false;
end

xDim=size(H,2);
zDim=size(H,1);

PMLInv=H'*pinv(R)*H;
PML=pinv(PMLInv);%The ML covariance matrix if R is known.
xEst=(PMLInv\H')*(R\z);%The ML estimate.

%If we should just use the ML estimate
if(xDim<=2)
    return;
end
%Otherwise, compute the shrinkage estimator.

%Equation 8
pStar=trace(PML)/max(eig(PML));

if(scaleUnknown~=false)    
    diff=z-H*xEst;
    %Equation 9
    sigma2=(diff'*pinv(R)*diff)/(zDim-xDim+2);
else
    sigma2=1;
end

%Equation 7, the shrinkage estimator.
xEst=max(0,1-sigma2*max(0,min((xDim-2),2*(pStar-2)))/(xEst'*PMLInv*xEst))*xEst;
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
