classdef InvGaussianD
%Functions to handle the scalar inverse Gaussian distribution. This
%distribution is used to describe the amount of time it takes a scalar
%Brownian motion model to reach a fixed value
%Implemented methods are: mean, var, PDF, CDF, momentGenFun, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(mu)
%%MEAN Obtain the mean of the inverse Gaussian distribution.
%
%INPUTS: mu The (positive) mean of the PDF. This is a scalar value.
%
%OUTPUTS: val The mean of the inverse Gaussian distribution.
%
%The inverse Gaussian distribution is parameterized by its mean and a shape
%parameter. Thus, this function just returns the parameter it is
%given.
%
%The inverse Gamma distribution is dicussed in [1]. 
%
%REFERENCES:
%[1] J. L. Folks and R. S. Chhikara, "The inverse Gaussian distribution and
%    its statistical application- a review," Journal of the Royal
%    Statistical Society. Series B. (Methodological), vol. 40, no. 3, pp.
%    263-289, 1978.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=mu;
end
   
function val=var(mu,lambda)
%%VAR   Obtain the variance of the inverse Gaussian distribution.
%
%INPUTS:    mu The mean of the PDF. This must be positive.
%       lambda The shape parameter of the PDF. This must be positive.
%
%OUTPUTS: val  The variance of the inverse Gaussian distribution.
%
%The variance of the inverse Gaussian distribution is given in [1].
%
%REFERENCES:
%[1] J. L. Folks and R. S. Chhikara, "The inverse Gaussian distribution and
%    its statistical application- a review," Journal of the Royal
%    Statistical Society. Series B. (Methodological), vol. 40, no. 3, pp.
%    263-289, 1978.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=mu^3/lambda;
end

function val=PDF(x,mu,lambda)
%%PDF Evaluate the inverse Gaussian probability distribution function (PDF)
%     at one or more desired points.
%
%INPUTS: x The point(s) at which the inverse Gaussian PDF is to be
%          evaluated.
%       mu The mean of the PDF. This must be positive.
%   lambda The shape parameter of the PDF. This must be positive.
%
%OUTPUTS: val The value(s) of the inverse Gaussian PDF.
%
%The PDF of the inverse Gaussian distribution is given in Equation 2 of
%[1].
%
%REFERENCES:
%[1] J. L. Folks and R. S. Chhikara, "The inverse Gaussian distribution and
%    its statistical application- a review," Journal of the Royal
%    Statistical Society. Series B. (Methodological), vol. 40, no. 3, pp.
%    263-289, 1978.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=sqrt(lambda./(2*pi*x.^3)).*exp(-lambda*(x-mu).^2./(2*mu^2*x));
    val(x<0)=0;
end

function val=CDF(x,mu,lambda)
%%CDF Evaluate the cumulative distribution function (CDF) of the inverse
%     Gaussian distribution at one or more desired points.
%
%INPUTS: x The point(s) at which the inverse Gaussian CDF is to be 
%          evaluated.
%       mu The mean of the PDF. This must be positive.
%   lambda The shape parameter of the PDF. This must be positive.
%
%OUTPUTS: val The value(s) of the inverse Gaussian CDF.
%
%Equation 7 in [1] expressed the CDF of the inverse Gaussian distribution
%in terms of the CDF of the standard normal distribution.
%
%EXAMPLE:
% mu=5;
% lambda=1;
% x=linspace(0,50,100);
% vals=InvGaussianD.CDF(x,mu,lambda);
% figure(1)
% clf
% plot(x,vals,'-k','linewidth',2)
%
%REFERENCES:
%[1] J. L. Folks and R. S. Chhikara, "The inverse Gaussian distribution and
%    its statistical application- a review," Journal of the Royal
%    Statistical Society. Series B. (Methodological), vol. 40, no. 3, pp.
%    263-289, 1978.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=GaussianD.CDF(sqrt(lambda./x).*(x/mu-1))+exp(2*lambda/mu).*GaussianD.CDF(-sqrt(lambda./x).*(x/mu+1));
end

function momentVal=momentGenFun(mu,lambda,numDerivs,t)
%%MOMENTGENFUN Evaluate the moment generating function (or one of its
%              derivatives) of the inverse Gaussian normal distribution.
%              Taking the ith derivative of the moment generating function
%              with respect to t and evaluating it at t=0 provides the
%              ith-order noncentral moment of the inverse Gaussian
%              distribution.
%
%INPUTS: mu The mean of the PDF. This must be positive.
%    lambda The shape parameter of the PDF. This must be positive.
% numDerivs The number of derivatives to take with respect to the argument
%           of the moment generating function. numDerivs>=0.
%         t The numPointsX1 or 1XnumPoints vector of points where the
%           moment generating function should be evaluated. If this
%           parameter is omitted or an empty matrix is passed, the default
%           value of 0 is used.
%
%OUTPUTS: momentVal A numPointsX1 vector of the values of the derivatives
%                   of the moment generating function given at the points
%                   in t or at t=0 if t is omitted.
%
%The moment generating function of a random variable is defined to be
%E(exp(t*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter. It can be shown that the moment
%Generating function of the inverse Gaussian distribution is
%f^0(t)=exp((lambda/mu)*(1-sqrt(1-2*mu^2*t/lambda)))
%This is related to the characteristic function, which is given in [1].
%This function evaluates the derivatives recursively using the chain rule
%and the fact that the nth derivative fo the argument fo the exponential
%term is easy to explicitly express. The characteristic function
%
%REFERENCES:
%[1] J. L. Folks and R. S. Chhikara, "The inverse Gaussian distribution and
%    its statistical application- a review," Journal of the Royal
%    Statistical Society. Series B. (Methodological), vol. 40, no. 3, pp.
%    263-289, 1978.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<4)
    t=0;
end

numT=length(t(:));
tVals=t;

momentVal=zeros(numT,1);
for curT=1:numT
    t=tVals(curT);
    
    %Derivatives of the argument of the exponent appear in the expression for
    %the derivative of the entire function.
    expArgDerivs=zeros(numDerivs+1,1);%Extra oen for zeroth-order derivative.
    expArgDerivs(0+1)=(lambda/mu)*(1-sqrt(1-(2*mu^2*t/lambda)));
    n=0:(numDerivs-1);
    expArgDerivs(2:end)=(mu.^(2*n+1).*doubleFactorial(2*n-1))./(lambda.^n.*(1-(2*t*mu^2)/lambda).^((2*n+1)/2));

    %The derivatives are constructed recursively using the chain rule, so we
    %have to keep track of all of the past ones.
    MGFDerivs=zeros(numDerivs+1,1);

    %The zeroth-order derivative is the function value itself.
    MGFDerivs(0+1)=exp(expArgDerivs(0+1));
    %Subsequent derivatives are recursively computed.

    %Precompute the binomial coefficients.
    binomTable=makeBinomTable(numDerivs-1);
    for curDeriv=1:numDerivs
        sumVal=0;
        for expDerivNum=0:(curDeriv-1)
            sumVal=sumVal+binomTable(curDeriv-1+1,expDerivNum+1)*MGFDerivs(expDerivNum+1)*expArgDerivs(curDeriv-expDerivNum+1);
        end
        MGFDerivs(curDeriv+1)=sumVal;
    end
    momentVal(curT)=MGFDerivs(end);
end
end

function val=rand(N,mu,lambda)
%%RAND Generate inverse Gaussian distributed random variables with the
%      given parameters.
%
%INPUTS:  N If N is a scalar, then rand returns an NXN matrix of random
%           variables. If N=[M,N1] is a two-element row vector, then rand
%           returns an MXN1 matrix of random variables.
%        mu The mean of the PDF. This must be positive.
%    lambda The shape parameter of the PDF. This must be positive.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated inverse Gaussian random variables.
%
%The random variables are generated using the algorithm given in [1].
%
%EXAMPLE:
%Here, we demonstrate that a histogram of many random samples will tend to
%agree with the PDF.
% mu=4;
% lambda=2;
% numSamp=100000;
% samp=InvGaussianD.rand([numSamp,1],mu,lambda);
% figure(1)
% clf
% hold on
% h=histogram(samp,'Normalization','pdf');
% %We will plot the PDF.
% sigma=sqrt(InvGaussianD.var(mu,lambda));
% numPoints=100;
% x=linspace(0,mu+4*sigma,numPoints);
% PDFVals=InvGaussianD.PDF(x,mu,lambda);
% plot(x,PDFVals,'-r','linewidth',2)
%
%REFERENCES:
%[1] J. R. Michael, W. R. Schucany, and R. W. Haas, "Generating random
%    variates using transformations with multiple roots," The American
%    Statistician, vol. 30, no. 2, pp. 88-90, May 1976.
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    v=randn(dims).^2;
    
    x1=mu*(1+(1/(2*lambda))*(mu*v-sqrt((4*mu*lambda+mu^2*v).*v)));
    
    u=rand(dims);
    
    val=x1;
    sel=(u>=mu./(mu+x1));
    val(sel)=mu^2./x1(sel);
end
    
end
end
