function [aEst,bEst]=binomialProbLikeBetaConjUpdate(xMeas,nMeas,aEst,bEst)
%%BINOMIALPROBLIKEBETACONJUPDATE When estimating the probability parameter
%           of the binomial distribution, the conjugate prior distribution
%           is a beta distribution. Given the prior shape parameters, this
%           function updates those parameters conditioned on the
%           measurements. The result is the parameters of the posterior
%           distribution, which is also a beta distribution.
%
%INPUTS:xMeas An NX1 or 1XN set of N independent measurements of the
%             binomial distribution.
%       nMeas The number of trials used in the binomial distribution for
%             all of the samples in xMeas.
%        aEst The shape parameter of the prior (beta) distribution that is
%             the exponent of x; aEst>0.
%        bEst The shape parameter of the prior distribution that is the
%             exponent of 1-x; bEst>0.
%
%OUTPUTS: aEst The shape parameter of the posterior (beta) distribution
%              that is the exponent of x. 
%         bEst The shape parameter of the posterior distribution that is
%              the exponent of 1-x.
%
%A distribution that is conjugate prior to a particular likelihood function
%is such that after performing Bayes' rule to update using a measurement,
%one obtains the same type of distribution back again. This function just
%performs Bayes rule with the given measurements.
%
%The conjugate prior relation for the negative binomial distribution is
%given in [1]. This can be derived in a similar manner.
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

aEst=aEst+sum(xMeas);
bEst=bEst+numMeas*nMeas-sum(xMeas);
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
