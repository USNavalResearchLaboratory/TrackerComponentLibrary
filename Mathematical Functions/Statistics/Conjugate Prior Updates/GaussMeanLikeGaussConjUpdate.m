function [muEst,SigmaEstInv]=GaussMeanLikeGaussConjUpdate(xMeas,SigmaMeasInv,muEst,SigmaEstInv)
%%GAUSSMEANLIKEGAUSSCONJUPDATE  When estimating the mean of a multivariate
%             normal distribution where the measurements have a known
%             covariance matrix, the conjugate prior is also a multivariate
%             normal distribution. Given the mean and covariance matrix of
%             the conjugate prior, as well as the covariance matrix of the
%             normal distribution whose mean is being esimated, this
%             function updates the mean and covariance matrix conditioned
%             on the measurements. These are the parameters of the
%             posterior distribution, which is also multivariate normal.
%
%INPUTS: xMeas An xDimXN set of N independent measurements of the
%              multivariate normal distribution with unknown mean and
%              inverse covariance matrix SigmaMeasInv.
% SigmaMeasInv The xDimXxDim inverse of the covariance matrix of the
%              likelihood function.
%        muEst The xDimX1 mean of the prior distribution, which is a
%              multivariate Gaussian distribution.
%  SigmaEstInv The xDimXxDim inverse of the covariance matrix of the prior
%              distribution, which is a multivariate Gaussian distribution.
%
%OUTPUTS: muEst The mean of the posterior distribution, which is a
%               multivariate normal distribution.
%   SigmaEstInv The inverse of the covariance matrix of the posterior
%               distribution.
%
%Note that while both SigmaMeasInv and SigmaEstInv can be singular, the
%sum, SigmaMeasInv+SigmaEstInv should be positive definite.
%
%A distribution that is conjugate prior to a particular likelihood function
%is such that after performing Bayes' rule to update using a measurement,
%one obtains the same type of distribution back again. This function just
%performs Bayes rule with the given measurements.
%
%The relations for the conjugate prior update are often homework
%assignments in statistics classes. A derivation is given in Chapter 2.3.3
%of [1].
%
%REFERENCES:
%[1] C. M. Bishop, Pattern Recognition and Machine Learning. Cambridge,
%    United Kingdom: Springer, 2007.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=size(xMeas,2);

SigmaEstInvPrev=SigmaEstInv;

SigmaEstInv=SigmaEstInv+numMeas*SigmaMeasInv;
muEst=SigmaEstInv\(SigmaEstInvPrev*muEst+SigmaMeasInv*sum(xMeas,2));

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
