classdef ChiSquareD
%%CHISQUARED Functions to handle the central or non-central chi-squared
%            distributions. These distributions arise when summing the
%            squares of independent normal random variables. The result is
%            central if the mean of the normal random variables was zero.
%Implemented methods are: mean, var, PDF, CDF, invCDF (only for the central
%                         chi squared distribution), rand, (only for the
%                         central chi squared distribution) entropy
%
%DEPENDENCIES: GammaD.m
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(nu,lambda)
%%MEAN Obtain the mean of chi squared probability distribution.
%
%INPUTS: nu The number of degrees of freedom of the chi-squared
%           distribution. Note that nu>0.
%    lambda The non-centrality parameter of the distribution. In the
%           central chi-squared distribution, this is zero.
%
%OUTPUTS: val The mean of the chi squared distribution under consideration.
%
%The noncentral chi-squared distribution is descirbed in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(lambda))
        lambda=0;
    end

    val=nu+lambda;
end

function val=var(nu,lambda)
%%VAR Obtain the variance of chi squared probability distribution for the
%     given number of degrees of freedom.
%
%INPUTS: nu The number of degrees of freedom of the chi-square
%           distribution. Note that nu>0.
%    lambda The non-centrality parameter of the distribution. In the
%           central chi-squared distribution, this is zero.
%
%OUTPUTS: val The variance of the chi squared distribution under
%             consideration.
%
%The noncentral chi-squared distribution is descirbed in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(lambda))
        lambda=0;
    end

    val=2*nu+4*lambda;
end
    
function val=PDF(x,nu,lambda)
%%PDF Evaluate the chi squared probability distribution function (PDF) at
%     the desired points.
%
%INPUTS: x The point(s) at which the chi-squared PDF is to be evaluated.
%          Note that x>=0 for nonzero probabilities.
%       nu The number of degrees of freedom of the chi-square distribution.
%          Note that nu>0.
%   lambda The non-centrality parameter of the distribution. In the central
%          chi-squared distribution, this is zero.
%
%OUTPUTS: val The value(s) of the chi squared PDF evaluated at x.
%
%The noncentral chi-squared distribution is descirbed in [1].
%
%When lambda is very small but nonzero, the PDF function could run into
%problems due to the logarithm of lambda. The logarithms were used instead
%of a simpler form to lessen the numerical instability. If an invalid value
%for the noncentral distribution is obtained, it is assumed that it is due
%to using an extremely small nonzero value of lambda, and the result is
%approximated using lambda=0.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(lambda))
        lambda=0;
    end

    if(lambda==0)
        val=(1/(2^(nu/2)*gamma(nu/2))).*x.^((nu-2)/2).*exp(-x/2);
    else
        val=0.5*exp(-(x+lambda)/2+(nu/4-0.5)*(log(x)-log(lambda))).*besseli(nu/2-1,sqrt(lambda*x));
        
        if(~isfinite(val))
            val=(1/(2^(nu/2)*gamma(nu/2))).*x.^((nu-2)/2).*exp(-x/2);
        end
    end
    
    val(x<0)=0;
end

function prob=CDF(x,nu,lambda)
%%CDF Evaluate the cumulative distribution function (CDF) of the central or
%     noncentral chi-squared distribution at desired points.
%
%INPUTS: x The point(s) at which the chi-squared CDF is to be evaluated.
%          Note that x>0 for nonzero values.
%       nu The number of degrees of freedom of the chi-square distribution.
%          Note that nu>0. When lambda is not zero, the CDF function is
%          numerically stabler when nu is an integer.
%   lambda The non-centrality parameter of the distribution. In the central
%          chi-squared distribution, this is zero.
%
%OUTPUTS: prob The value(s) of the CDF of the chi-squared distribution with
%              nu degrees of freedom and noncentrality parameter lambda
%              evaluated at x.
%
%The CDF of the central chi squared function with nu degrees of freedom
%evaluated at x is just a special case of the incomplete gamma function, as
%described at [1]. The incomplete gamma function is built into Matlab
%without the use of any toolboxes. When the noncentrality parameter is
%nonzero, the CDF is evaluated in terms of the equivalent noncentral gamma
%distribution. The noncentral gamma algorithm is documented further in
%GammaD.m.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%Updated to use the noncentral gamma distribution's function
%December 2014 David A. Karnick, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(lambda))
        lambda=0;
    end

    if(lambda==0)
        prob=gammainc(x/2,nu/2);
    else
        prob=GammaD.CDF(x,nu/2,2,lambda);
    end
    prob(x<=0)=0;
end
    
function val=invCDF(prob,nu)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the central chi squared distribution. This is only implemented
%        for the central chi-squared distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%          nu The number of degrees of freedom of the chi-squared
%             distribution. Note that nu>0.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The CDF of the chi squared function with nu degrees of freedom evaluated
%at x is just a special case of the incomplete gamma function, as described
%on Mathworld at [1]. The inverse incomplete gamma function is built into
%Matlab without the use of any toolboxes and thus is called here.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=2*gammaincinv(prob,nu/2);
end

function vals=rand(N,nu,lambda)
%%RAND Generate chi-squared distributed random variables with the given
%      parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%       nu The number of degrees of freedom of the chi-squared
%          distribution. Note that nu>0. If lambda is not zero, nu >=1.
%   lambda The non-centrality parameter of the distribution. In the central
%          chi-squared distribution, this is zero. If this parameter is
%          omitted or an empty matrix is passed, then 0 is used.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated chi-squared random variables.
%
%The algorithm for the central chi squared distribution is an
%implementation of the inverse transform algorithm of Chapter 5.1 of [1].
%When the noncentral distribution is used, the random variables are
%generated by summing the squares of normally distributed random variables.
%
%A non-central chi squared distribution is the sum of a chi-squared
%distribution with nu-1 degrees of freedom plus the squared of a normal
%distribution with mean sqrt(lambda) and standard deviation 1. This is how
%the noncentral distribution is generated and why it is required that
%nu>=1.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(lambda))
        lambda=0;
    end

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    if(lambda==0)
        U=rand(dims);
        vals=ChiSquareD.invCDF(U,nu);
    else
        U=rand(dims);
        vals=ChiSquareD.invCDF(U,nu-1)+(sqrt(lambda)+randn(N)).^2;
    end
end

function entropyVal=entropy(nu)
%%ENTROPY Obtain the differential entropy of the central chi squared
%         distribution given in nats. The differential entropy of a
%         continuous distribution is entropy=-int_x p(x)*log(p(x)) dx where
%         the integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: nu The number of degrees of freedom of the chi-squared
%           distribution. Note that nu>0 
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
    
    entropyVal=nu/2+log(2)+gammaln(nu/2)+(1-nu/2)*psi(nu/2);
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
