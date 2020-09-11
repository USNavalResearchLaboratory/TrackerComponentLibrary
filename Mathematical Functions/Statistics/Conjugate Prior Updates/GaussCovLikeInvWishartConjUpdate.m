function [PsiEst,nuEst]=GaussCovLikeInvWishartConjUpdate(xMeas,muMeas,PsiEst,nuEst)
%%GAUSSCOVLIKEINVWISHARTCONJUPDATE When estimating the covariance matrix of
%             a multivariate normal distribution where the measurements
%             have a known mean, the conjugate prior is an inverse Wishart
%             distribution. Given the precision matrix and number of
%             degrees of freedom of the inverse Wishart distribution, as
%             well as the mean of the normal distribution whose covariance
%             matrix is being estimated, this function updates the
%             precision matrix and number of degrees of freedom conditioned
%             on the measurements. These are the parameters of the
%             posterior distribution, which is also inverse Wishart.
%
%INPUTS: xMeas An xDimXN set of N independent measurements of the
%              multivariate normal distribution with unknown mean and
%              inverse covariance matrix SigmaMeasInv.
%       muMeas The xDimX1 mean of the likelihood function.
%       PsiEst The xDimXxDim positive-definite, symmetric precision matrix
%              of the prior (inverse Wishart) distribution.
%        nuEst The number of degrees of freedom of the prior distribution.
%              Note that nu>=xDim.
%
%OUTPUTS: PsiEst The precision matrix of the posterior distribution, which
%                is an inverse Wishart distribution.
%          nuEst The number of degrees of freedom of the posterior
%                distribution.
%
%A distribution that is conjugate prior to a particular likelihood function
%is such that after performing Bayes' rule to update using a measurement,
%one obtains the same type of distribution back again. This function just
%performs Bayes rule with the given measurements.
%
%The conjugate prior update is given in Chapter 7.7.2 of [1].
%
%REFERENCES:
%[1] T. W. Anderson, An Introduction to Multivariate Statistical Analysis,
%    3rd ed. Hoboken, NJ: Wiley-Interscience, 2003.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=size(xMeas,2);

nuEst=nuEst+numMeas;

diff=bsxfun(@minus,xMeas,muMeas);
PsiEst=PsiEst+diff*diff';

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
