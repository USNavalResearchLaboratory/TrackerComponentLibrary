function SStar=covMatEstShrinkage(x,centerData,algorithm)
%%COVMATESTSHRINKAGE When estimating a covariance matrix from a small
%           number of samples, it can occur that the resultng estimate is
%           poorly conditioned. This function estimates the covariance
%           matrix from samples using one of two shrinkage estimators that
%           are mean to avoid poor conditioning and that can produce more
%           accurate (in terms of the squared Frobenius norm) estimates
%           than one can obtain from the sample covariance matrix.
%
%INPUTS: x A pXn matrix of n samples of a p-dimensional distribution.
% centerData If the samples in x are zero mean or the mean has been
%          subtracted from them, then centerData can be false. Otherwise,
%          centerData should be true to center the data (subtract the
%          sample mean). The algorithmic choices assume that the true mean
%          has been subtracted, so subtracting the sample mean might affect
%          the performance. The default if omitted or an empty matrix is
%          passd is true.
% algorithm An optional parameter selecting the algorithm to use. Possible
%          values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            Ledoit-Wolf algorithm of [1].
%          1 Use the empirical Bayesian estimator (the Hanff estimator),
%            which is also described in [1].
%
%OUTPUTS: SStar The pXp covariance matrix estimate.
%
%Note that the algorithms are not unbiased. However, the Ledoit-Wolf
%algorithm is shown in [1] to be asymptotically unbiased (as n gets large
%keeping p/n fixed).
%
%EXAMPLE:
%Here, we look at the performance when considering random diagonal
%covariance matrices whose elements are drawn from the log-normal
%distribution. We look at the percentage relative improvement in average
%loss (PRIAL) criterion, which is essentially the relative improvement in
%the squared Frobenius norm of the error of the sample covariance matrix
%versus the Ledoit-Wolf covariance matrix esitmate.
% numRuns=10000;
% p=20;
% n=40;
% %Parameters for the log-normal distribution from which the covariance
% %matrix is drawn.
% mu=1;
% Sigma=1;
% centerData=false;
%
% SSigmaFro2=0;
% SStarSigmaFro2=0;
% for curRun=1:numRuns
%     x=LogNormalD.rand(p,mu,Sigma);
%     X=diag(x);%The current covariance matrix.
%     SX=cholSemiDef(X,'lower');
%     x=SX*randn(p,n);%Zero mean.
%     
%     S=(1/n)*(x*x');%Sample covariance matrix.
%     SStar=covMatEstShrinkage(x,centerData);
% 
%     SSigmaFro2=SSigmaFro2+norm(S-Sigma,'fro')^2;
%     SStarSigmaFro2=SStarSigmaFro2+norm(SStar-Sigma,'fro')^2;
% end
% SSigmaFro2=SSigmaFro2/numRuns;
% SStarSigmaFro2=SStarSigmaFro2/numRuns;
% PRIAL=(SSigmaFro2-SStarSigmaFro2)/SSigmaFro2
%One will typically see that the PRIAL is positive, around 0.08, indicating
%that the Ledoit-Wolf covariance matrix esitmate is 8% better than the
%sample covariance estimate. It was empirically noticed that setting
%centerData=true will further improve the performance on this example even
%though the samples generated are already zero-mean.
%
%REFERENCES:
%[1] O. Ledoit and M. Wolf "A Well-conditioned estimator for large-
%    dimensional covariance matrices," Journal of Multivariate Analysis,
%    vol. 88, no. 2, pp. 365-411, 2004.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(algorithm))
    algorithm=0;
end

if(nargin<2||isempty(centerData))
    centerData=true;
end

if(centerData)
    x=bsxfun(@minus,x,mean(x,2)); 
end

p=size(x,1);
n=size(x,2);
%The sample covariance matrix.
S=(1/n)*(x*x');

switch(algorithm)
    case 0%The Ledoit-Wolf covariance estimator.
        %Lemma 3.2
        m=trace(S)/p;

        %Lemma 3.3
        d2=norm(S-diag(m),'fro')^2;

        %Lemma 3.4
        b2Bar=0;
        for k=1:n
            b2Bar=b2Bar+norm(x(:,k)*x(:,k)'-S,'fro')^2;
        end
        b2Bar=(1/n^2)*b2Bar;
        b2=min(b2Bar,d2);

        %Lemma 3.5
        a2=d2-b2;

        %Equation 14
        SStar=(b2/d2)*diag(m)+(a2/d2)*S;
    case 1%The empirical Bayesian estimator (Haff estimator).
        mEB=det(S)^(1/p);
        if(mEB<=0||~isreal(mEB))
            mEB=trace(S)/p;
        end

        SStar=((p*n-2*n-2)/(p*n^2))*diag(mEB)+(n/(n+1))*S;
    otherwise
        error('Unknown algorithm specified.')
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
