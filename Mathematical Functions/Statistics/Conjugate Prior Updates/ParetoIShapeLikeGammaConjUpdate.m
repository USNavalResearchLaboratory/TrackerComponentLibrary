function [kEst,thetaEst]=ParetoIShapeLikeGammaConjUpdate(xMeas,xMinMeas,kEst,thetaEst)
%%PARETOISHAPELIKEGAMMACONJUPDATE When estimating the shape parameter of
%                  the Pareto type I distribution with a known scale
%                  parameter, the conjugate prior is a central (lambda=0)
%                  gamma distribution. Given the shape and scale parameters
%                  of the conjugate prior as well as the scale parameter of
%                  the measurement (Pareto) distribution, this function
%                  updates the parameters of the prior distribution
%                  conditioned on the measurements. The result is the
%                  parameters of the posterior distribution,
%                  which is also a gamma distribution.
%
%INPUTS:xMeas An NX1 or 1XN set of N independent measurements of the
%              uniform distribution.
%      xMinEst The scale parameter of the likelihood distribution (a Pareto
%              type I distribution). This marks the lowest point having
%              nonzero likelihood. This must be >0.
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
%The conjugate prior relation for the Pateto type I distribution is given
%in [1].
%
%REFERENCES:
%[1] B. C. Arnold and S. J. Press, "Bayesian estimation and prediction for
%    Pareto data," Journal of the American Statistical Association, vol.
%    84, no. 408, pp. 1079-1084, Dec. 1989.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numMeas=length(xMeas);

kEst=kEst+numMeas;
betaEst=1/thetaEst;
thetaEst=1/(betaEst+sum(log(xMeas/xMinMeas)));

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
