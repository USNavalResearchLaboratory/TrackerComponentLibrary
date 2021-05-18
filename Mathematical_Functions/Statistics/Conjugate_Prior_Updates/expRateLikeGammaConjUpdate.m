function [kEst,thetaEst]=expRateLikeGammaConjUpdate(xMeas,kEst,thetaEst)
%%EXPRATELIKEGAMMACONJUPDATE When estimating the rate parameter of the
%               exponential distribution, the conjugate prior is the
%               central (lambda=0) gamma distribution. Given the shape
%               parameter and scale parameter of the prior central gamma
%               distribution,  this function updates those parameters using
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
