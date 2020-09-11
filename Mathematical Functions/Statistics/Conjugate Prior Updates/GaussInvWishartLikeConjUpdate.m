function [muUpdate,kappaUpdate,nuUpdate,VUpdate]=GaussInvWishartLikeConjUpdate(xList,mu,kappa,V,nu)
%%GAUSSINVWISHARTCONJUPDATE When estimating the mean and covariance matrix
%          of a multivariate Gaussian distribution, the conjugate prior
%          is a Gaussian- inverse-Wishart distribution, also called a
%          normal-inverse Wishart distribution. The normal inverse Wishart
%          PDF of of the form N(x;mu,Sigma/kappa)*IW(Sigma;nu,V), where N is
%          the normal PDF, IW is the inverse Wishart PDF, x and Sigma are
%          random variables and mu, kappa, nu, and V are parameters being
%          update dby this function given n samples of x.
%
%INPUTS: xList The xDimXn set of n measurements.
%           mu The xDimX1 mean of the  normal part of the prior
%              distribution.
%        kappa The positive, scalar scale factor of the normal part of the
%              distribution.
%            V The xDimXxDim precision matrix of the inverse Wishart
%              distribution.
%           nu The scalar number of degrees of freedom of the inverse
%              Wishart distribution; nu>d+1. 
%
%OUTPUTS: muUpdate, kappaUpdate, nuUpdate, VUpdate The updated parameters
%              of the Gaussian-inverse Wishart distribution.
%
%This conjugate prior relationship can play a role in extended target
%tracking algorithms.
%
%A distribution that is conjugate prior to a particular likelihood function
%is such that after performing Bayes' rule to update using a measurement,
%one obtains the same type of distribution back again. This function just
%performs Bayes rule with the given measurements.
%
%The conjugate prior relationship is derived in Chapter 3.6 of [1].
%
%REFERENCES:
%[1] A. Gelman, J. B. Carlin, S. H. S., D. B. Bunson, A. Vehtari, and D. B.
%    Rubin, Bayesian Data Analysis, 3rd ed. Boca Raton: CRC Press, 2013.
%
%May 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(xList,2);

[xBar,X]=calcMixtureMoments(xList);
X=n*X;

muUpdate=kappa/(kappa+n)*mu+n/(kappa+n)*xBar;
kappaUpdate=kappa+n;
nuUpdate=nu+n;
VUpdate=V+X+(kappa*n/(kappa+n))*(xBar-mu)*(xBar-mu)';

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
