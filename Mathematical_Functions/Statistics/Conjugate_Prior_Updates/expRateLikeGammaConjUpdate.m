function [kEst,thetaEst]=expRateLikeGammaConjUpdate(xMeas,kEst,thetaEst)
%%EXPRATELIKEGAMMACONJUPDATE When estimating the rate parameter of the
%               exponential distribution, the conjugate prior is the
%               central (lambda=0) gamma distribution. Given the shape
%               parameter and scale parameter of the prior central gamma
%               distribution, this function updates those parameters using
%               conditioned on the measurements. The result is the
%               parameters of the posterior distribution, which is also
%               cental gamma.
%
%INPUTS: xMeas An NX1 or 1XN set of N independent measurements of the
%              exponential distribution.
%         kEst The shape parameter of the prior (central gamma)
%              distribution.
%     thetaEst The scale parameter of the prior distribution.
%
%OUTPUTS: kEst The shape parameter of the posterior (central gamma)
%              distribution.
%     thetaEst The scale parameter of the posterior distribution.
%
%A distribution that is conjugate prior to a particular likelihood function
%is such that after performing Bayes' rule to update using a measurement,
%one obtains the same type of distribution back again. This function just
%performs Bayes rule with the given measurements.
%
%The conjugate prior relation for the exponential distribution is given in
%[1].
%
%EXAMPLE:
%The function needs a prior distribution to work. Here, we show how a prior
%distribution based on assumed minimum and maximum values for lambda
%suffices. The normalized estimation error squared (NEES) of the mean and
%variance comig from the estimates produced by this function (for the gamma
%distribution) are generally consistent (near 1).
% numMCRuns=1e4;
% minVal=10;
% maxVal=1000;
% %We will put a prior distribution on lambda. The mean will be set midway
% %between the minimum and maximum values and the variance will be
% %(maxVal-minVal)^2.
% %For the gamma distribution:
% %mean=k*theta
% %variance=k*theta^2
% %So, we are solving
% %k*theta=(minVal+maxVal)/2
% %k*theta^2=(maxVal-minVal)^2
% %So,
% kInit=(maxVal+minVal)^2/(4*(maxVal-minVal)^2);
% thetaInit=2*(maxVal-minVal)^2/(maxVal+minVal);
% 
% numMeas=10;
% NEES1=0;
% NEES2=0;
% for curMCRun=1:numMCRuns
%     %MUST have a prior. It will not work without a prior.
%     k=kInit;
%     theta=thetaInit;
%     for curMeas=1:numMeas
%         x1=ExponentialD.rand(1,minVal);
%         [k,theta]=expRateLikeGammaConjUpdate(x1,k,theta);
%     end
%     estMean=GammaD.mean(k,theta);
%     estVar=GammaD.var(k,theta);
%     NEES1=NEES1+(estMean-minVal)^2/estVar;
% 
%     k=kInit;
%     theta=thetaInit;
%     for curMeas=1:numMeas
%         x2=ExponentialD.rand(1,maxVal);
%         [k,theta]=expRateLikeGammaConjUpdate(x2,k,theta);
%     end
%     estMean=GammaD.mean(k,theta);
%     estVar=GammaD.var(k,theta);
%     NEES2=NEES2+(estMean-maxVal)^2/estVar;
% end
% NEES1=NEES1/numMCRuns
% NEES2=NEES2/numMCRuns
%
%REFERENCES:
%[1] D. Fink, "A compendium of conjugate priors," Montana State University,
%    Department of Biology, Environmental Statistics Group, Tech. Rep., May
%    1997. [Online].
%    Available: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.157.5540
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=length(xMeas);

kEst=kEst+numMeas;

%The inverse of the scale parameter=rate parameter of the gamma
%distribution.
beta=1/thetaEst;
thetaEst=1./(beta+sum(xMeas));

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
