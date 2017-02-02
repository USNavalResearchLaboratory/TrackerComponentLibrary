function [kEst,betaEst]=gammaScaleLikeInvGammaConjUpdate(xMeas,kMeas,kEst,betaEst)
%%GAMMASCALELIKECONJUPDATE When estimating the scale parameter of a central
%                   (lambda=0) gamma distribution having a known shape
%                   parameter, the conjugate prior distribution is an
%                   inverse gamma distribution. Given the shape parameter
%                   of the measurement distribution and the prior shape and
%                   rate parameter of the inverse gamma distribution, this
%                   function updates the prior parameters conditioned on
%                   the measurements. The result is the parameters of the
%                   posterior distribution, which is also an inverse gamma
%                   distribution.
%
%INPUTS: xMeas An NX1 or 1XN set of N independent measurements of the
%              central gamma distribution with shape parameter kMeas.
%        kMeas The shape parameter of the likelihood function.
%        kEst  The shape parameter of the prior distribution, which is an
%              inverse gamma distribution.
%      betaEst The rate parameter of the prior distribution.
%
%OUTPUTS:kEst The shape parameter of the posterior distribution, which is
%             an inverse gamma distribution.
%     betaEst The rate parameter of the posterior distribution.
%
%A distribution that is conjugate prior to a particular likelihood function
%is such that after performing Bayes' rule to update using a measurement,
%one obtains the same type of distribution back again. This function just
%performs Bayes rule with the given measurements.
%
%The conjugate prior relation for the gamma distribution with known shape
%parameter is given in [1].
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

kEst=kEst+numMeas*kMeas;
betaEst=betaEst+sum(xMeas);
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
