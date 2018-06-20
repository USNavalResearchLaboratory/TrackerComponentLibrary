function [muEst,deltaEst]=gammaShapeLikeGammaConjIConjUpdate(xMeas,thetaMeas,muEst,deltaEst)
%%GAMMASHAPELIKEGAMMACONJICONJUPDATE When estimating the shape parameter of
%               a gamma distribution having a known scale parameter, the
%               conjugate prior distribution is known as the gamma
%               conjugate type I distribution. Given the scale parameter of
%               the measurement distribution and the two parameters of the
%               prior distribution, this function updates the parameters of
%               the prior distribution conditioned on the measurements. The
%               result is the parameters of the posterior distribution,
%               which is also a gamma conjugate type I distribution.
%
%INPUTS: xMeas An NX1 or 1XN set of N independent measurements of the
%              central gamma distribution with scale parameter bMeas.
%    thetaMeas The known scale parameter of the likelihood of the
%              measurements; thetaMeas>0.
%        muEst The parameter of the prior distribution (gamma conjugate
%              type I)that is raised to delta*x; mu>0.
%     deltaEst The parameter of the prior distribution that is the exponent
%              of gamma(x); delta>0.
%
%OUTPUTS:muEst The parameter of the posterior distribution (gamma conjugate
%              type I)that is raised to delta*x; mu>0.
%     deltaEst The parameter of the posterior distribution that is the
%              exponent of gamma(x); delta>0.
%
%A distribution that is conjugate prior to a particular likelihood function
%is such that after performing Bayes' rule to update using a measurement,
%one obtains the same type of distribution back again. This function just
%performs Bayes rule with the given measurements.
%
%The conjugate prior relation for the gamma distribution with known scale
%parameter is given in Theorem 1 in [1].
%
%REFERENCES:
%[1] E. Damsleth, "Conjugate classes for gamma distributions," Scandinavian
%    Journal of Statistics, vol. 2, no. 2, pp. 80-84, 1975.  
%   
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=length(xMeas);

muS=muEst*thetaMeas;
muS=exp(log(muS)*(deltaEst/(deltaEst+numMeas))+log(geometricMean(xMeas))*(numMeas/(deltaEst+numMeas)));
muEst=muS/thetaMeas;

deltaEst=deltaEst+numMeas;

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
