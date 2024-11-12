classdef ChiSquareD
%%CHISQUARED Functions to handle the central or non-central chi-squared
%            distributions. These distributions arise when summing the
%            squares of independent normal random variables. The result is
%            central if the mean of the normal random variables was zero.
%
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand, entropy (only
%                         for the central chi squared distribution),
%                         nthMoment
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(nu,lambda,alpha)
%%MEAN Obtain the mean of chi squared probability distribution.
%
%INPUTS: nu The number of degrees of freedom of the chi-squared
%           distribution. Note that nu>0.
%    lambda The non-centrality parameter of the distribution. In the
%           central chi-squared distribution, this is zero. If this
%           parameter is omitted or an empty matrix is passed, then 0 is
%           used.
%     alpha A scale factor that appears in the definition of the noncentral
%           chi squared distribution in [2] and represents a slightly
%           nonstandard variant of the chi-squared distribution. The
%           default if omitted or an empty matrix is passed is 1, which is
%           the value in the more standard chi squared distribution in [1].
%
%OUTPUTS: val The mean of the chi squared distribution under consideration.
%
%The noncentral chi-squared distribution is described in [1]. The version
%with the alpha term is used in [2].
%
%EXAMPLE:
%This shows the sample mean of the noncentral chi squared distribution
%and the mean returned by this function are close.
% numSamp=1e6;
% nu=8;
% lambda=16;
% samples=ChiSquareD.rand([1,numSamp],nu,lambda);
% sampleMean=mean(samples)
% trueMean=ChiSquareD.mean(nu,lambda)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%[2] M. S. Holla, "On a noncentral chi-square distribution in the analysis
%    of weapon systems effectiveness," Metrika: International Journal for
%    Theoretical and Applied Statistics, vol. 15, pp. 9-14, Dec. 1970.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(alpha))
        alpha=1;
    end

    if(nargin<2||isempty(lambda))
        lambda=0;
    end

    val=(nu+lambda)/alpha;
end

function val=var(nu,lambda,alpha)
%%VAR Obtain the variance of noncentral chi squared probability
%     distribution.
%
%INPUTS: nu The number of degrees of freedom of the chi-square
%           distribution. Note that nu>0.
%    lambda The non-centrality parameter of the distribution. In the
%           central chi-squared distribution, this is zero. If this
%           parameter is omitted or an empty matrix is passed, then 0 is
%           used.
%     alpha A scale factor that appears in the definition of the noncentral
%           chi squared distribution in [2] and represents a slightly
%           nonstandard variant of the chi-squared distribution. The
%           default if omitted or an empty matrix is passed is 1, which is
%           the value in the more standard chi squared distribution in [1].
%
%OUTPUTS: val The variance of the chi squared distribution under
%             consideration.
%
%The noncentral chi-squared distribution is described in [1]. The version
%with the alpha term is used in [2].
%
%EXAMPLE:
%This shows the sample variance of the noncentral chi squared distribution
%and the variance returned by this function are close.
% numSamp=1e6;
% nu=8;
% lambda=16;
% samples=ChiSquareD.rand([1,numSamp],nu,lambda);
% sampleVar=var(samples)
% trueVar=ChiSquareD.var(nu,lambda)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%[2] M. S. Holla, "On a noncentral chi-square distribution in the analysis
%    of weapon systems effectiveness," Metrika: International Journal for
%    Theoretical and Applied Statistics, vol. 15, pp. 9-14, Dec. 1970.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(alpha))
        alpha=1;
    end

    if(nargin<2||isempty(lambda))
        lambda=0;
    end

    val=2*(nu+2*lambda)/alpha^2;
end

function val=nthMoment(n,nu,lambda,alpha,nMaxHyper)
%%NTHMOMENT Compute nth Obtain the varnoncentral moment of the noncentral 
%     chi squared probability distribution. n need not be an integer.
%
%INPUTS: n The order of the moment. n>=0.
%       nu The number of degrees of freedom of the chi-square
%          distribution. Note that nu>0.
%   lambda The non-centrality parameter of the distribution. In the central
%          chi-squared distribution, this is zero. If this parameter is
%          omitted or an empty matrix is passed, then 0 is used.
%    alpha A scale factor that appears in the definition of the noncentral
%          chi squared distribution in [2] and represents a slightly
%          nonstandard variant of the chi-squared distribution. The
%          default if omitted or an empty matrix is passed is 1, which is
%          the value in the more standard chi squared distribution in [1].
% nMaxHyper This parameter can usually be omitted. If n is not an integer
%          or n>10, then an expression utilizing hypergeometric1F1 is
%          used. This corresponds to the nMax input in that function. If
%          omitted or an empty matrix is passed, the default in the
%          hypergeometric1F1 function is used.
%
%OUTPUTS: val The nth noncentral moments of the chi squared distribution
%             under consideration.
%
%The noncentral chi-squared distribution is described in [1]. The version
%with the alpha term is used in [2]. The moments were found analytically.
%For n being an integer between 0 and 10, an explicit expression is used.
%For n>=11 or being a non-integer, an solution using hypergeometric1F1 is
%used.
%
%EXAMPLE:
%In this example, one can see the number of digits of agreement between the
%nth sample moment for the noncentral chi squared distribution and the
%value returned by this function.
% numSamp=1e6;
% n=1.2345;
% nu=5;
% lambda=16;
% samples=ChiSquareD.rand([1,numSamp],nu,lambda);
% sampleMoment=mean(samples.^n)
% trueMoment=ChiSquareD.nthMoment(n,nu,lambda)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%[2] M. S. Holla, "On a noncentral chi-square distribution in the analysis
%    of weapon systems effectiveness," Metrika: International Journal for
%    Theoretical and Applied Statistics, vol. 15, pp. 9-14, Dec. 1970.
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<5||isempty(nMaxHyper))
    nMaxHyper=[];
end

if(nargin<4||isempty(alpha))
    alpha=1;
end

if(nargin<3||isempty(lambda))
    lambda=0;
end

if(n<0)
    error('n must be >=0')
end

if(imag(n)~=0)
    error('n must be real.')
end

switch(n)
    case 0
        val=1;
    case 1
        val=(nu+lambda)/alpha;
    case 2
        val=(nu^2+2*nu*(1+lambda)+lambda*(4+lambda))/alpha^2;
    case 3
        val=(nu*(2+nu)*(4+nu)+3*(2+nu)*(4+nu)*lambda+3*(4+nu)*lambda^2+lambda^3)/alpha^3;
    case 4
        val=(nu*(2+nu)*(4+nu)*(6+nu)+4*(2+nu)*(4+nu)*(6+nu)*lambda+6*(4+nu)*(6+nu)*lambda^2+4*(6+nu)*lambda^3+lambda^4)/alpha^4;
    case 5
        val=1/alpha^5*(nu*(2+nu)*(4+nu)*(6+nu)*(8+nu)+5*(2+nu)*(4+nu)*(6+nu)*(8+nu)*lambda+10*(4+nu)*(6+nu)*(8+nu)*lambda^2+10*(6+nu)*(8+nu)*lambda^3+5*(8+nu)*lambda^4+lambda^5);
    case 6
        val=1/alpha^6*(nu*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)+6*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)*lambda+15*(4+nu)*(6+nu)*(8+nu)*(10+nu)*lambda^2+20*(6+nu)*(8+nu)*(10+nu)*lambda^3+15*(8+nu)*(10+nu)*lambda^4+6*(10+nu)*lambda^5+lambda^6);
    case 7
        val=1/alpha^7*(nu*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)+7*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*lambda+21*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*lambda^2+35*(6+nu)*(8+nu)*(10+nu)*(12+nu)*lambda^3+35*(8+nu)*(10+nu)*(12+nu)*lambda^4+21*(10+nu)*(12+nu)*lambda^5+7*(12+nu)*lambda^6+lambda^7);
    case 8
        val=1/alpha^8*(nu*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)+8*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*lambda+28*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*lambda^2+56*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*lambda^3+70*(8+nu)*(10+nu)*(12+nu)*(14+nu)*lambda^4+56*(10+nu)*(12+nu)*(14+nu)*lambda^5+28*(12+nu)*(14+nu)*lambda^6+8*(14+nu)*lambda^7+lambda^8);
    case 9
        val=1/alpha^9*(nu*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)+9*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)*lambda+36*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)*lambda^2+84*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)*lambda^3+126*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)*lambda^4+126*(10+nu)*(12+nu)*(14+nu)*(16+nu)*lambda^5+84*(12+nu)*(14+nu)*(16+nu)*lambda^6+36*(14+nu)*(16+nu)*lambda^7+9*(16+nu)*lambda^8+lambda^9);
    case 10
        val=1/alpha^10*(nu*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)*(18+nu)+10*(2+nu)*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)*(18+nu)*lambda+45*(4+nu)*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)*(18+nu)*lambda^2+120*(6+nu)*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)*(18+nu)*lambda^3+210*(8+nu)*(10+nu)*(12+nu)*(14+nu)*(16+nu)*(18+nu)*lambda^4+252*(10+nu)*(12+nu)*(14+nu)*(16+nu)*(18+nu)*lambda^5+210*(12+nu)*(14+nu)*(16+nu)*(18+nu)*lambda^6+120*(14+nu)*(16+nu)*(18+nu)*lambda^7+45*(16+nu)*(18+nu)*lambda^8+10*(18+nu)*lambda^9+lambda^10);
    otherwise
        val=2^n*gamma(nu/2+n)*alpha^(-n)*hypergeometric1F1(-n,nu/2,-lambda/2,nMaxHyper)/gamma(nu/2);
end
end

function val=PDF(x,nu,lambda,alpha,AbsTol,RelTol)
%%PDF Evaluate the chi squared probability distribution function (PDF) at
%     the desired points.
%
%INPUTS: x The point(s) at which the chi-squared PDF is to be evaluated.
%          Note that x>=0 for nonzero probabilities.
%       nu The number of degrees of freedom of the chi-square distribution.
%          Note that nu>0.
%   lambda The non-centrality parameter of the distribution. In the central
%          chi-squared distribution, this is zero. If this parameter is
%          omitted or an empty matrix is passed, then 0 is used.
%    alpha A scale factor that appears in the definition of the noncentral
%          chi squared distribution in [2] and represents a slightly
%          nonstandard variant of the chi-squared distribution. The
%          default if omitted or an empty matrix is passed is 1, which is
%          the value in the more standard chi squared distribution in [1].
% AbsTol, RelTol If lambda is not zero, then the ChiSquaredSumD.PDF is used
%          and these correspond to the same-named inputs of that function.
%          The defaults if omitted or empty matrices are passed are 1e-11
%          and 1e-7.
%
%OUTPUTS: val The value(s) of the chi squared PDF evaluated at x.
%
%The noncentral chi-squared distribution is described in [1].
%
%When lambda is zero, then the evaluate of the PDF is generally
%straightforward. When lambda is not zero, then the PDF expression in [1]
%involving the besseli function can easily have overflow problems. Thus,
%when lambda is not zero, the integration method of ChiSquaredSumD.PDF is
%used.
%
%EXAMPLE:
%This plots a histogram of random samples of the noncentral chi squared
%distribution along with the PDF. One can see that the results agree well.
% numSamp=1e5;
% alpha=2.123;
% nu=8;
% lambda=1000;
% samples=ChiSquareD.rand([1,numSamp],nu,lambda,alpha);
% numPts=1000;
% x=linspace(0,600,numPts);
% y=ChiSquareD.PDF(x,nu,lambda,alpha);
% figure(1)
% clf
% hold on
% histogram(samples,'normalization','pdf')
% plot(x,y,'linewidth',2);
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%[2] M. S. Holla, "On a noncentral chi-square distribution in the analysis
%    of weapon systems effectiveness," Metrika: International Journal for
%    Theoretical and Applied Statistics, vol. 15, pp. 9-14, Dec. 1970.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<6||isempty(RelTol))
        RelTol=1e-7;
    end

    if(nargin<6||isempty(AbsTol))
        AbsTol=1e-11;
    end

    if(nargin<4||isempty(alpha))
        alpha=1;
    end

    if(nargin<3||isempty(lambda))
        lambda=0;
    end

    x=alpha*x;

    if(lambda==0)
        val=(1/(2^(nu/2)*gamma(nu/2))).*x.^((nu-2)/2).*exp(-x/2);
    else
        val=ChiSquaredSumD.PDF(x,1,nu,lambda,0,0,AbsTol,RelTol);
        %The above is more stable than using something like.
        %val=0.5*exp(-(x+lambda)/2+(nu/4-0.5)*(log(x)-log(lambda))).*besseli(nu/2-1,sqrt(lambda*x));
    end
    
    val=alpha*val;
    val(x<0)=0;
end

function prob=CDF(x,nu,lambda,alpha,AbsTol,RelTol)
%%CDF Evaluate the cumulative distribution function (CDF) of the central or
%     noncentral chi-squared distribution at desired points.
%
%INPUTS: x The point(s) at which the chi-squared CDF is to be evaluated.
%          Note that x>0 for nonzero values.
%       nu The number of degrees of freedom of the chi-square distribution.
%          Note that nu>0. When lambda is not zero, the CDF function is
%          numerically stabler when nu is an integer.
%   lambda The non-centrality parameter of the distribution. In the central
%          chi-squared distribution, this is zero. If this parameter is
%          omitted or an empty matrix is passed, then 0 is used.
%    alpha A scale factor that appears in the definition of the noncentral
%          chi squared distribution in [2] and represents a slightly
%          nonstandard variant of the chi-squared distribution. The
%          default if omitted or an empty matrix is passed is 1, which is
%          the value in the more standard chi squared distribution in [1].
% AbsTol, RelTol If lambda is not zero, then the ChiSquaredSumD.CDF is used
%          and these correspond to the same-named inputs of that function.
%          The defaults if omitted or empty matrices are passed are 1e-11
%          and 1e-7.
%
%OUTPUTS: prob The value(s) of the CDF of the chi-squared distribution with
%              nu degrees of freedom and noncentrality parameter lambda
%              evaluated at x.
%
%The CDF of the central chi squared function with nu degrees of freedom
%evaluated at x is just a special case of the incomplete gamma function, as
%described at [1]. The incomplete gamma function is built into Matlab
%without the use of any toolboxes. When the noncentrality parameter is
%nonzero, the CDF could be evaluated in terms of the equivalent noncentral
%gamma distribution. However, it was found that ChiSquaredSumD.CDF is much
%more numericallty stable for large lambda.
%
%EXAMPLE:
%In this example, we generate a bunch of random samples and display the
%empirical CDF. We also display the theoretical CDF. One can see that they
%agree well.
% numSamp=1e5;
% alpha=2.1;
% nu=8;
% lambda=1600;
% samples=ChiSquareD.rand([1,numSamp],nu,lambda,alpha);
% numPts=2000;
% x=linspace(0,1400,numPts);
% CDFEmp=EmpiricalD.CDF(x,samples);
% CDF=ChiSquareD.CDF(x,nu,lambda,alpha);
% figure(1)
% clf
% hold on
% plot(x,CDFEmp,'linewidth',4)
% plot(x,CDF,'linewidth',2)
% legend('Empirical CDF','Theoretical CDF','location','northwest')
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%[2] M. S. Holla, "On a noncentral chi-square distribution in the analysis
%    of weapon systems effectiveness," Metrika: International Journal for
%    Theoretical and Applied Statistics, vol. 15, pp. 9-14, Dec. 1970.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%Updated to use the noncentral gamma distribution's function
%December 2014 David A. Karnick, Naval Research Laboratory, Washington D.C.
%Updated to use the ChiSquaredSumD.CDF function instead of the gamma
%function.
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<6||isempty(RelTol))
        RelTol=1e-11;
    end

    if(nargin<6||isempty(AbsTol))
        AbsTol=1e-7;
    end

    if(nargin<4||isempty(alpha))
        alpha=1;
    end

    if(nargin<3||isempty(lambda))
        lambda=0;
    end

    x=x*alpha;

    if(lambda==0)
        prob=gammainc(x/2,nu/2);
    else
        %This is less subject to finite precision issues than GammaD.CDF.
        prob=ChiSquaredSumD.CDF(x,1,nu,lambda,0,0,AbsTol,RelTol);
        %prob=GammaD.CDF(x,nu/2,2,lambda);
    end
    prob(x<=0)=0;
end
    
function val=invCDF(prob,nu,lambda,alpha,otherParams)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the chi squared distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%         argument of the CDF is desired.
%      nu The number of degrees of freedom of the chi-squared distribution.
%          Note that nu>0.
%   lambda The non-centrality parameter of the distribution. In the central
%          chi-squared distribution, this is zero. If this parameter is
%          omitted or an empty matrix is passed, then 0 is used.
%    alpha A scale factor that appears in the definition of the noncentral
%          chi squared distribution in [2] and represents a slightly
%          nonstandard variant of the chi-squared distribution. The
%          default if omitted or an empty matrix is passed is 1, which is
%          the value in the more standard chi squared distribution in [1].
% otherParams If lamdba is not zero, then this optinal parameter is a
%          structure that can have members named 'numStdInt', 'BrentVals',
%          and 'intVals', which correspond to the same-named inputs of
%          ChiSquaredSumD.invCDF. See the comments to ChiSquaredSumD.invCDF
%          for more information. If omitted or an empty matrix is provided,
%          the defaults in ChiSquaredSumD.invCDF are used.s
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The CDF of the chi squared function with nu degrees of freedom and
%lambda=0 evaluated at x is just a special case of the incomplete gamma
%function, as described on Mathworld at [1]. The inverse incomplete gamma
%function is built into Matlab without the use of any toolboxes and thus is
%called here. For the case of lambda being nonzero, ChiSquaredSumD.invCDF
%is called and the input otherParams is used, if provided.
%
%EXAMPLE:
%Here, for a given distribution, we get the probability at a particular
%point and then show that when using invCDF, one can get the same point
%back. Also, we start with a given probability and use invCDF to get a
%point.
% alpha=2.123;
% nu=8;
% lambda=1000;
% xTrue=450;
% prob=ChiSquareD.CDF(xTrue,nu,lambda,alpha);
% xBack=ChiSquareD.invCDF(prob,nu,lambda,alpha);
% RelxErr=(xBack-xTrue)/xTrue
% 
% probTrue=0.999;
% x=ChiSquareD.invCDF(probTrue,nu,lambda,alpha);
% probBack=ChiSquareD.CDF(x,nu,lambda,alpha);
% RelProbErr=(probBack-probTrue)/probTrue
%
%REFERENCES:
%[1] Weisstein, Eric W. "Noncentral Chi-Squared Distribution." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/NoncentralChi-SquaredDistribution.html
%[2] M. S. Holla, "On a noncentral chi-square distribution in the analysis
%    of weapon systems effectiveness," Metrika: International Journal for
%    Theoretical and Applied Statistics, vol. 15, pp. 9-14, Dec. 1970.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(alpha))
        alpha=1;
    end

    if(nargin<3||isempty(lambda))
        lambda=0;
    end

    if(lambda==0)
        val=2*gammaincinv(prob,nu/2);
    else
        numStdInt=[];
        BrentVals=[];
        intVals=[];

        if(nargin>4&&~isempty(otherParams))
            if(isfield(otherParams,'numStdInt'))
                numStdInt=otherParams.numStdInt;
            end
            if(isfield(otherParams,'BrentVals'))
                BrentVals=otherParams.BrentVals;
            end
            if(isfield(otherParams,'intVals'))
                intVals=otherParams.intVals;
            end
        end

        val=ChiSquaredSumD.invCDF(prob,1,nu,lambda,0,0,numStdInt,BrentVals,intVals);
    end
    val=val/alpha;
end

function vals=rand(N,nu,lambda,alpha,algorithm)
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
%    alpha A scale factor that appears in the definition of the noncentral
%          chi squared distribution in [2] and represents a slightly
%          nonstandard variant of the chi-squared distribution. The
%          default if omitted or an empty matrix is passed is 1, which is
%          the value in the more standard chi squared distribution in [1].
% algorithm An optional parameter that changes how the chi squared random
%          variables are generated. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use the
%            inverse transform algorithm of Chapter 5.1 of [1], as
%            described below.
%          1 This can only be chosen if nu is an integer. Sum up the
%            squares of nu normal random variables to get the appropriate
%            chi squared random variable from the relation between normal
%            and chi squared random variables.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated chi-squared random variables.
%
%In algorithm 0, the algorithm for the central chi squared distribution is
%an implementation of the inverse transform algorithm of Chapter 5.1 of
%[1]. A non-central chi squared distribution is the sum of a chi-squared
%distribution with nu-1 degrees of freedom plus the squared of a normal
%distribution with mean sqrt(lambda) and standard deviation 1. This is how
%the noncentral distribution is generated and why it is required that
%nu>=1.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%[2] M. S. Holla, "On a noncentral chi-square distribution in the analysis
%    of weapon systems effectiveness," Metrika: International Journal for
%    Theoretical and Applied Statistics, vol. 15, pp. 9-14, Dec. 1970.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<5||isempty(algorithm))
    algorithm=0;
end

if(nargin<4||isempty(alpha))
    alpha=1;
end

if(nargin<3||isempty(lambda))
    lambda=0;
end

if(isscalar(N))
    dims=[N, N];
else
    dims=N;
end

if(algorithm==0)
    if(lambda==0)
        U=rand(dims);
        vals=ChiSquareD.invCDF(U,nu,0,alpha);
    else
        U=rand(dims);
        vals=(ChiSquareD.invCDF(U,nu-1)+(sqrt(lambda)+randn(N)).^2)/alpha;
    end
elseif(algorithm==1)
    if(fix(nu)~=nu)
        error('Algorithm 1 only works when nu is an integer.')
    end

    vals=zeros(dims);
    numEls=prod(dims);
    mu=sqrt(lambda/nu);
    vals(:)=sum((mu+randn(nu,numEls)).^2,1)/alpha;
else
    error('Unknown algorithm specified.')
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
