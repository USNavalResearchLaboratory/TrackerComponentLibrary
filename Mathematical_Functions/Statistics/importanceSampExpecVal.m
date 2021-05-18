function [mu,P]=importanceSampExpecVal(h,f,g,gSamp,numTrials)
%%IMPORTANCESAMPEXPECVAL Approximate the expected value of h(x) with
%               respect to the probability density function (PDF) f(x)
%               using Monte Carlo integration with importance sampling.
%               Importance sampling is useful if a direct sampling of f is
%               difficult and/ or if the variance of h(x) is large.
%               Importance sampling samples a distribution g instead of f
%               and converges faster than standard Monte-Carlo integration
%               if the variance of h(x)f(x)/g(x) is small.
%
%INPUTS: h A handle to a function (possibly multivariate) over which the
%          expected value is to be estimated. The output should be a column
%          vector.
%        f A function handle for the (possibly multivariate) PDF with
%          respect to which the integral over h is taken.
%        g A function handle for the PDF that is easy to sample. This must
%          have support of at least the intersection of the support of f
%          and g.
%    gSamp A function handle to draw samples of the distribution h.
%          hSamp() should return a single random sample of h as a
%          numDimX1 vector.
%numTrials The number of Monte Carlo trials to perform for the numerical
%          integration. If this parameter is omitted or an empty matrix is
%          passed, numTrials=1000 is used.
%
%OUTPUTS: mu The approximate expected value of h(x) with respect to the PDF
%            f(x).
%          P An unbiased estimate of the sample covariance matrix of mu.
%            This will be an empty matrix if numTrials=1, because no
%            unbiased estimate is available. 
%
%Importance sampling is described in Chapter 8.6 of [1]. The basic idea of
%Monte Carlo integration is discussed in Chapter 3.2 of [1]. The matrix P
%comes from the fact that each sample has the covariance matrix associated
%with h(x)f(x)/g(x) and the covariance matrix of the average is 1/numTrials
%times that (from Monte Carlo integration). Thus, P is computed by taking
%the unbiased sample average of the covariance matrix of h(x)f(x)/g(x) and
%dividing it by numTrials.
%
%EXAMPLE:
%Suppose we want to evaluate the integral from 0 to 1 of exp(x^2). In this
%case f(x)=1 for x between 0 and 1. We choose g(x)=(1/(exp(1)-1))*exp(x)
%between 0 and 1 and 0 elsewhere. This means that
%h(x)=(exp(1)-1)*exp(x^2-x) and leads to a CDF of
%prob(x)=(exp(x)-1)/(exp(1)-1). The inverse of this CDF is
%x=log(1+(exp(1)-1)*prob). Thus we can sample g(x) using an inverse
%transform sampling algorithm (See Chapter 5.1 of [1]) to get gSamp.
% numTrials=10000;
% h=@(x)(exp(1)-1)*exp(x*(x-1));
% f=@(x)1;
% g=@(x)(1/(exp(1)-1))*exp(x);
% gSamp=@()log(1+(exp(1)-1)*rand(1));
% [muIS,PIS]=importanceSampExpecVal(h,f,g,gSamp,numTrials)
% %In comparison, a Monte Carlo approximation is 
% U=rand(1,numTrials);x=exp(U.^2);
% muMC=mean(x)
% PMC=std(x)^2/numTrials
%One should see that the variance PIS of the importance sampling estimate
%is less than that of the standard Monte Carlo estimate PMC and that muIS
%is usually closer to the true value of
%mu=1.4626517459071816088040485868570 than muMC.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(numTrials))
    numTrials=1000;
end

if(nargout<2)
    %There is no need to save all of the samples if only the mean is
    %desired.
    mu=0;
    for curSamp=1:numTrials
        x=gSamp();
        mu=mu+h(x)*f(x)/g(x);
    end
    mu=mu/numTrials;
else
   %Save the samples so a sample mean and covariance matrix can be found. 
   x=gSamp();
   mu1=h(x)*f(x)/g(x);

   muDim=length(mu1);
   muSamp=zeros(muDim,numTrials,1);
   muSamp(:,1)=mu1;
   for curSamp=2:numTrials
        x=gSamp();
        muSamp(:,curSamp)=h(x)*f(x)/g(x);
   end
   mu=mean(muSamp,2);
   
   if(numTrials==1)
       P=[];%No unbiased covariance estimate available.
   else
       P=zeros(muDim,muDim);
       for curSamp=1:numTrials
           diff=muSamp(:,curSamp)-mu;
           P=P+diff*diff';
       end
       P=P/((numTrials-1)*numTrials);
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
