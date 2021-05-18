classdef ChiD
%%CHID Functions to handle the chi distribution. 
%These distributions arise when taking the square root of the sum of the
%squares of independent normal random variables.
%Implemented methods are: mean, var, moments, PDF, CDF, (supports
%                         noncentral distributions) rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(nu)
%%MEAN Obtain the mean of central chi probability distribution.
%
%INPUTS: nu The number of degrees of freedom of the central chi
%           distribution. Note that nu>0.
%
%OUTPUTS: val The mean of the central chi distribution.
%
%The mean of the central chi distribution is given in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Chi Distribution." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/ChiDistribution.html
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    val=sqrt(2)*exp(gammaln((nu+1)/2)-gammaln(nu/2));
end

function val=var(nu)
%%VAR Obtain the variance of the central chi probability distribution.
%
%INPUTS: nu The number of degrees of freedom of the chi distribution. Note
%           that nu>0.
%
%OUTPUTS: val  The variance of the central chi distribution.
%
%The variance of the central chi distribution is given in [1]. 
%
%REFERENCES:
%[1] Weisstein, Eric W. "Chi Distribution." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/ChiDistribution.html
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=2*exp(gammaln(1+nu/2)-gammaln(nu/2))-ChiD.mean(nu)^2;
end

function val=moments(m,nu)
%%MOMENTS Obtain noncentral moments of the chi distribution.
%
%INPUTS: m The order of the moment m>=0. The zeroth order moment is just 1.
%       nu The number of degrees of freedom of the chi distribution. Note
%          that nu>0.
%
%OUTPUTS: val  The mth-degree noncentral moments of the chi distribution.
%
%Moments of the chi distribution are given in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Chi Distribution." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/ChiDistribution.html
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=2^(m/2)*exp(gammaln((1/2)*(nu+m))-gammaln(nu/2));
end

function val=PDF(x,nu)
%%PDF Evaluate the central chi probability distribution function (PDF).
%
%INPUTS: x The point(s) at which the central chi PDF is to be evaluated.
%          Note that x>=0 for a nonzero probability.
%       nu The number of degrees of freedom of the chi distribution. Note
%          that nu>0. Numerical precision problems can arise for very large
%          values of nu.
%
%OUTPUTS: val The value(s) of the central chi PDF evaluated at x.
%
%The PDF is given in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Chi Distribution." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/ChiDistribution.html
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(x));
    sel=(x>=0);
    
    y=x(sel);
    val(sel)=2^(1-nu/2)*y.^(nu-1).*exp(-y.^2/2)/gamma(nu/2);
end

function prob=CDF(x,nu)
%%CDF Evaluate the cumulative distribution function (CDF) of the central
%     chi distribution at desired points.
%
%INPUTS: x The point(s) at which the chi CDF is to be evaluated.
%       nu The number of degrees of freedom of the chi-square
%          distribution. Note that nu>0.
%
%OUTPUTS: prob The value(s) of the CDF of the chi distribution with nu
%              degrees of freedom evaluated at x.
%
%The CDF of the central chi distribution is given in [1]. 
%
%REFERENCES:
%[1] Weisstein, Eric W. "Chi Distribution." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/ChiDistribution.html
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    prob=gammainc(x.^2/2,nu/2);
end

function vals=rand(N,nu,lambda)
%%RAND Generate chi distributed random variables with the given parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%       nu The number of degrees of freedom of the chi distribution. Note
%          that nu>0. If lambda is not zero, nu >=1.
%   lambda The non-centrality parameter of the distribution. In the central
%          chi distribution, this is zero. If this parameter is omitted or
%          an empty matrix is passed, then 0 is used.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated chi-squared random variables.
%
%A chi distributed random variable is just the square root of a chi-squared
%random variable. Thus, this function takes the square root of
%ChiSquareD.rand(N,nu,lambda).
%
%EXAMPLE:
%Here, we demonstrate that a histogram of many random samples will tend to
%agree with the PDF of the central distribution.
% nu=4.4;
% lambda=0;
% numSamp=[10000,1];
% samp=ChiD.rand(numSamp,nu);
% figure(1)
% clf
% hold on
% h=histogram(samp,'Normalization','pdf');
% %We will plot the PDF.
% numPoints=100;
% x=linspace(0,20,numPoints);
% PDFVals=ChiD.PDF(x,nu);
% plot(x,PDFVals,'-r','linewidth',2)
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(lambda))
        lambda=0;
    end

    if(isscalar(N))
        dims=[N N];
    else
        dims=N;
    end

    vals=sqrt(ChiSquareD.rand(dims,nu,lambda));
end

function entropyVal=entropy(nu)
%%ENTROPY Obtain the differential entropy of the central chi distribution
%         given in nats. The differential entropy of a continuous
%         distribution is entropy=-int_x p(x)*log(p(x))  dx where the
%         integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: nu The number of degrees of freedom of the chi distribution. Note
%           that nu>0. 
%
%OUTPUTS: entropyVal The value of the differential entropy in nats.
%
%Differential entropy is defined in Chapter 8 of [1].
%
%REFERENCES:
%[1] T. M. Cover and J. A. Thomas, Elements of Information Theory, 2nd ed.
%    Hoboken, NJ: Wiley-Interscience, 2006.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    entropyVal=gammaln(nu/2)+(1/2)*(nu-log(2)-(nu-1)*psi(nu/2)); 
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
