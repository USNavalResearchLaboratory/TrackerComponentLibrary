classdef GammaConjugateID
%%GAMMACONJUGATEID Functions to handle the gamma conjugate type I
%          distribution. This distribution is introduced in [1] as the
%          conjugate prior to the central gamma distribution with a known
%          scale parameter when the shape parameter is being estimated.
%Implemented methods are: mean, PDF, normConst
%
%REFERENCES:
%[1] E. Damsleth, "Conjugate classes for gamma distributions," Scandinavian
%    Journal of Statistics, vol. 2, no. 2, pp. 80-84, 1975.  
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
 
methods(Static)
function val=mean(mu,delta,normConst,varargin)
%%MEAN Evaluate the mean of the gamma conjugate type I distribution for a
%      given set of parameters.
%
%INPUTS:    mu The parameter of the distribution that is raised to delta*x;
%              mu>0.
%        delta The parameter of the distribution that is the exponent of
%              gamma(x); delta>0. If delta is large (a few hundred), then
%              finite precision errors are likely to cause problems with
%              the numeric integration.
%    normConst Optionally the normalizing constant of the distribution can 
%              be provided. If omitted or an empty matrix is passed,
%              the function GamconID.normConst is used to obtain it. This
%              is offered as an option, because the computation of the
%              normalizing constant requires numerical integration and can
%              be slow.
%     varargin Any parameters that should be passed to Matlab's integral
%              function to control the precision. This is generally a
%              series of text strings followed by values.
%
%OUTPUTS: val The mean of the gamma conjugate type I distribution.
%
%The distribution is described in [1]. However, there is no closed-form
%solution for the mean. Thus, the mean is computed using numerical
%integration via Matlab's integral function.
%
%REFERENCES:
%[1] E. Damsleth, "Conjugate classes for gamma distributions," Scandinavian
%    Journal of Statistics, vol. 2, no. 2, pp. 80-84, 1975.  
%   
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(normConst))
    normConst=GammaConjugateID.normConst(mu,delta);
end

%PDFScaled is the PDF, but with substitution so that the bounds of the
%integral can be from 0 to 1 instead of from 0 to Inf.
val=integral(@(y)fun2IntScaled(y,delta,mu,normConst),0,1,varargin{:});
function vals=fun2IntScaled(y,delta,mu,normConst)
    %This is equivalent to vals=normConst*gamma(y./(1-y)).^(-delta).*mu.^(delta*y./(1-y)).*(y./(1-y).^3);
    vals=exp(log(normConst)-delta*gammaln(y./(1-y))+(delta*y./(1-y))*log(mu)+log((y./(1-y).^3)));
end
end

function val=PDF(x,mu,delta,normConst)
%%PDF Evaluate the probability density function (PDF) of the gamma
%     conjugate type I distribution at one or more desired points.
%
%INPUTS:     x The point(s) at which the Weibull PDF is to be evaluated.
%           mu The parameter of the distribution that is raised to delta*x;
%              mu>0.
%        delta The parameter of the distribution that is the exponent of
%              gamma(x); delta>0. If delta is large (a few hundred), then
%              finite precision errors are likely to cause problems with
%              the numeric integration if normConst needs to be computed.
%    normConst Optionally the normalizing constant of the distribution can 
%              be provided. If omitted or an empty matrix is passed,
%              the function GamconID.normConst is used to obtain it. This
%              is offered as an option, because the computation of the
%              normalizing constant requires numerical integration and can
%              be slow.
%
%OUTPUTS: val The value(s) of the gamma conjugate type I PDF with the given
%             parameters.
%
%The PDF of the gamma conjugate type I distribution is given in [1].
%
%REFERENCES:
%[1] E. Damsleth, "Conjugate classes for gamma distributions," Scandinavian
%    Journal of Statistics, vol. 2, no. 2, pp. 80-84, 1975.  
%   
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4||isempty(normConst))
    normConst=GammaConjugateID.normConst(mu,delta);
end

%The following is equivalent to val=normConst*gamma(x).^(-delta).*mu.^(delta*x);
val=exp(log(normConst)-delta*gammaln(x)+log(mu)*(delta*x));
val(x<0)=0;
end

function val=normConst(mu,delta,varargin)
%%NORMCONST Compute the normalization cconstant for the gamma conjugate
%           type I distribution. This multiplied by the unnormalized PDF
%           results in a normalized PDF.
%
%INPUTS:    mu The parameter of the distribution that is raised to delta*x;
%              mu>0.
%        delta The parameter of the distribution that is the exponent of
%              gamma(x); delta>0. If delta is large (a few hundred), then
%              finite precision errors are likely to cause problems with
%              the numeric integration.
%     varargin Any parameters that should be passed to Matlab's integral
%              function to control the precision. This is generally a
%              series of text strings followed by values.
%
%OUTPUTS: val The normalization constant.
%
%As there is no analytic method for computing the normalization constant,
%numerical integration over the PDF, which is given in [1], is performed
%using Matlab's integral function.
%
%REFERENCES:
%[1] E. Damsleth, "Conjugate classes for gamma distributions," Scandinavian
%    Journal of Statistics, vol. 2, no. 2, pp. 80-84, 1975.  
%   
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

%PDFScaled is the PDF, but with substitution so that the bounds of the
%integral can be from 0 to 1 instead of from 0 to Inf.
val=1/integral(@(y)PDFScaled(y,delta,mu),0,1,varargin{:});
function vals=PDFScaled(y,delta,mu)
    %The following is equivalent to gamma(y./(1-y)).^(-delta).*mu.^(delta*y./(1-y)).*(1./(1-y).^2);
    vals=exp(-delta*gammaln(y./(1-y))+log(mu)*(delta*y./(1-y))-2*log(1-y));
end
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
