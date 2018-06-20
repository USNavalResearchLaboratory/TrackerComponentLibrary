classdef KD
%%KD Functions to handle the scalar K-distribution. The K-distribution is
%    often used to model the observed radar cross section (RCS) of radar
%    returns due to clutter, as mentioned in [1]. Multiple definitions of
%    _the_ K distribution exist depending on how many parameters there are.
%    This function implements the distributions for one, two, and three
%    parameters. These are often called the K0, K and K' distributions (or
%    all just the K distribution) as in [2].
%Implemented methods are: mean, var, PDF, moments, rand
%
%REFERENCES:
%[1] N. J. Redding, "Estimating the parameters of the K distribution in
%    the intensity domain," DSTO Electronics and Surveillance Research
%    Laboratory, Salisbury, South Australia, Australia, Tech. Rep. DTSOTR-
%    0839, Jul. 1999.
%[2] M. C. Teich and P. Diament, "Multiply stochastic representations for
%    K distributions and their Poisson transforms," Journal of the Optical
%    Society of America A, vol. 6, no. 1, pp. 80-91, Jan. 1989.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(mu)
%%MEAN Obtain the mean of the K distribution.
%
%INPUTS:   mu  The (positive) mean of the PDF. This is a scalar value.
%
%OUTPUTS: val  The mean of the K distribution.
%
%The K distribution is parameterized by its mean and two other
%parameters. Thus, this function just returns the parameter it is
%given.
%
%Moments of the K distribution are shown in Table 1 in [1].
%
%REFERENCES:
%[1] M. C. Teich and P. Diament, "Multiply stochastic representations for
%    K distributions and their Poisson transforms," Journal of the Optical
%    Society of America A, vol. 6, no. 1, pp. 80-91, Jan. 1989.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=mu;
end
    
function val=var(mu,alpha,beta)
%%VAR Obtain the variance of the K distribution.
%
%INPUTS:    mu The (positive) mean of the PDF. This is a scalar value.
%  alpha, beta The degrees of freedom parameters for the two gamma
%              distributions that make up the K distribution. The K
%              distribution is obtained from the integral 
%       int_0^inf GammaD.PDF(x,beta,y/beta)*GammaD.PDF(y,alpha,mu/alpha) dy
%              These are the parameters of the teo Gamma distributions. mu
%              is the mean of the innermost one and alpha is the number of
%              degrees of freedom.
%
%OUTPUTS: val  The variance of the K distribution.
%
%Moments of the K distribution are shown in Table 1 in [1].
%
%REFERENCES:
%[1] M. C. Teich and P. Diament, "Multiply stochastic representations for
%    K distributions and their Poisson transforms," Journal of the Optical
%    Society of America A, vol. 6, no. 1, pp. 80-91, Jan. 1989.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(mu^2*(1+alpha+beta))/(alpha*beta);
end

function val=PDF(x,mu,alpha,beta)
%%PDF Obtain theprobability density function (PDF) of the K distribution
%
%INPUTS:     x The points at which the PDF should be evaluated. This can be
%              a scalar or a matrix.
%           mu The (positive) mean of the PDF. This is a scalar value.
%  alpha, beta The degrees of freedom parameters for the two gamma
%              distributions that make up the K distribution. The K
%              distribution is obtained from the integral 
%       int_0^inf GammaD.PDF(x,beta,y/beta)*GammaD.PDF(y,alpha,mu/alpha) dy
%              These are the parameters of the teo Gamma distributions. mu
%              is the mean of the innermost one and alpha is the number of
%              degrees of freedom.
%
%OUTPUTS: val  The PDF of the K distribution evaluated at the points in x.
%
%PDFs of the K distribution are shown in Table 1 in [1]. THe K0
%distribution has alpha=1 and beta=1. The standard K distiribution has just
%beta=1. Here, we allow for alpha and beta to be positive, values,
%representing the general K' distribution in [1].
%
%REFERENCES:
%[1] M. C. Teich and P. Diament, "Multiply stochastic representations for
%    K distributions and their Poisson transforms," Journal of the Optical
%    Society of America A, vol. 6, no. 1, pp. 80-91, Jan. 1989.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3)
    %The K0 distribution in [2], which can be seen to be
    %equivalent to the K distribution in [2] if alpha=1.
        alpha=1;
    end

    if(nargin<4)
    %The K distribution in [2], which is equivalent to the K'
    %distribution with beta=1.
        beta=1;
    end

    val=2*(alpha*beta/mu)*(1/(gamma(alpha)*gamma(beta)))*(alpha*beta*x/mu).^((alpha+beta)/2-1).*besselk(alpha-beta,2*sqrt(alpha*beta*x/mu));
    %Negative points have zero probability as this is a one-sided
    %distribution.
    val(x<=0)=0;
end
    
function val=moments(m,mu,alpha,beta)
%%MOMENTS Obtain noncentral moments of the K distribution.
%
%INPUTS:     m The order of the noncentral moment m>=0. The zeroth order
%              moment is just 1.
%           mu The (positive) mean of the PDF. This is a scalar value.
%  alpha, beta The degrees of freedom parameters for the two gamma
%              distributions that make up the K distribution. The K
%              distribution is obtained from the integral 
%       int_0^inf GammaD.PDF(x,beta,y/beta)*GammaD.PDF(y,alpha,mu/alpha) dy
%              These are the parameters of the teo Gamma distributions. mu
%              is the mean of the innermost one and alpha is the number of
%              degrees of freedom.
%
%OUTPUTS: val  The mth-degree noncentral moments of the K distribution.
%
%Moments of the K distribution are shown in Table 1 in [1].
%
%REFERENCES:
%[1] M. C. Teich and P. Diament, "Multiply stochastic representations for
%    K distributions and their Poisson transforms," Journal of the Optical
%    Society of America A, vol. 6, no. 1, pp. 80-91, Jan. 1989.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(mu/(alpha*beta))^m*risingFactorial(alpha,m)*risingFactorial(beta,m);
end

function vals=rand(N,mu,alpha,beta)
%%RAND       Generate K distributed random variables with the given
%            parameters.
%
%INPUTS:    N  If N is a scalar, then rand returns an NXN matrix of random
%              variables. If N=[M,N1] is a two-element row vector, then
%              rand returns an MXN1 matrix of random variables.
%           mu The (positive) mean of the PDF. This is a scalar value.
%  alpha, beta The degrees of freedom parameters for the two gamma
%              distributions that make up the K distribution. The K
%              distribution is obtained from the integral 
%       int_0^inf GammaD.PDF(x,beta,y/beta)*GammaD.PDF(y,alpha,mu/alpha) dy
%              These are the parameters of the teo Gamma distributions. mu
%              is the mean of the innermost one and alpha is the number of
%              degrees of freedom.
%
%OUTPUTS:   vals   A matrix whose dimensions are determined by N of the
%                  generated K random variables.
%
%The K distribution is discussed in [1]. We generating samples by sampling
%one gamma distribution and them sampling another conditioned on the value
%obtained from the first. This function does not generate correlated K
%distributed random variables.
%
%EXAMPLE:
%Here, we demonstrate that a histogram of many random samples will tend to
%agree with the PDF.
% mu=10;
% alpha=4;
% beta=2;
% numSamp=[10000,1];
% samp=KD.rand(numSamp,mu,alpha,beta);
% figure(1)
% clf
% hold on
% h=histogram(samp,'Normalization','pdf');
% %We will plot the PDF.
% sigma=sqrt(KD.var(mu,alpha,beta));
% numPoints=100;
% x=linspace(0,mu+4*sigma,numPoints);
% PDFVals=KD.PDF(x,mu,alpha,beta);
% plot(x,PDFVals,'-r','linewidth',2)
%
%REFERENCES:
%[1] M. C. Teich and P. Diament, "Multiply stochastic representations for
%    K distributions and their Poisson transforms," Journal of the Optical
%    Society of America A, vol. 6, no. 1, pp. 80-91, Jan. 1989.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    %Gamma distribution scale parameter
    theta=mu/alpha;
    k=alpha;

    %The first distribution generates mean values that go into the
    %second distribution.
    vals=GammaD.rand(dims,k,theta);

    numVals=numel(vals);

    k=beta;
    for curVal=1:numVals
        theta=vals(curVal)/beta;

        vals(curVal)=GammaD.rand(1,k,theta);
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
