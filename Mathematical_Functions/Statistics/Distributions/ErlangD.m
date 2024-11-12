classdef ErlangD
%%ERLANGD Functions to handle the Erlang distribution.
%Implemented methods are: mean, var, nthMoment, PDF, CDF, invCDF,
%                         momentGenFun, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(k,lambda)
%%MEAN Obtain the mean of the Erlang distribution.
%
%INPUTS: k The shape parameter of the Erlang distribution (k>=1).
%   lambda The rate parameter of the Erlang distribution (lambda>0).
%
%OUTPUTS: val The mean of the exponential distribution.
%
%The PDF of the Erlang distribution is given in [1] and the mean can be
%analytically found by integrating x times the PDF.
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% k=3;
% lambda=1/4;
% meanVal=ErlangD.mean(k,lambda)
% numSamp=1e6;
% sampMeanVal=mean(ErlangD.rand([numSamp,1],k,lambda))
%One will find both values are about 12.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Erlang Distribution." From MathWorld--A Wolfram
%    Web Resource. https://mathworld.wolfram.com/ErlangDistribution.html
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=k/lambda;
end

function val=var(k,lambda)
%%VAR Obtain the variance of the Erlang distribution.
%
%INPUTS: k The shape parameter of the Erlang distribution (k>=1).
%   lambda The rate parameter of the Erlang distribution (lambda>0).
%
%OUTPUTS: val The variance of the Erlang distribution.
%
%The PDF of the Erlang distribution is given in [1] and the variance can be
%found analytically.
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% k=3;
% lambda=1/4;
% varVal=ErlangD.var(k,lambda)
% numSamp=1e6;
% sampVarVal=var(ErlangD.rand([numSamp,1],k,lambda))
%One will find both values are about 48.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Erlang Distribution." From MathWorld--A Wolfram
%    Web Resource. https://mathworld.wolfram.com/ErlangDistribution.html
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=k/lambda^2;
end

function val=nthMoment(n,k,lambda)
%%NTHMOMENT Obtain the nth noncentral moment of the Erlang distribution.
%
%INPUTS: n The moment number n>=0.
%        k The shape parameter of the Erlang distribution (k>=1).
%   lambda The rate parameter of the Erlang distribution (lambda>0).
%
%OUTPUTS: val The nth noncentral moment of the Erlang distribution.
%
%The PDF of the Erlang distribution is given in [1] and the nth noncentral
%moment can be found analytically by integrating x^n times the PDF.
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% n=5;
% k=3;
% lambda=1/4;
% varVal=ErlangD.nthMoment(n,k,lambda)
% numSamp=1e6;
% sampVarVal=mean(ErlangD.rand([numSamp,1],k,lambda).^n)
%One will find both values are relatively close in value.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Erlang Distribution." From MathWorld--A Wolfram
%    Web Resource. https://mathworld.wolfram.com/ErlangDistribution.html
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=gamma(k+n)/(lambda^n*gamma(k));
end

function val=PDF(x,k,lambda)
%%PDF Evaluate the PDF of the Erlang probability distribution function at
%     one or more desired points.
%
%INPUTS: x The point(s) at which the Erlang distribution is to be 
%          evaluated.
%        k The shape parameter of the Erlang distribution (k>=1).
%   lambda The rate parameter of the Erlang distribution (lambda>0).
%
%OUTPUTS: val The value(s) of the Erlang PDF.
%
%The PDF of the Erlang distribution is given in Chapter 4 of [1] and in
%[2].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% k=3;
% lambda=1/5;
% numSamples=1e5;
% figure(1)
% clf
% histogram(ErlangD.rand([numSamples,1],k,lambda),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,60,numPoints);
% vals=ErlangD.PDF(x,k,lambda);
% plot(x,vals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%[2] Weisstein, Eric W. "Erlang Distribution." From MathWorld--A Wolfram
%    Web Resource. https://mathworld.wolfram.com/ErlangDistribution.html
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(lambda*(lambda*x).^(k-1)/gamma(k)).*exp(-lambda*x);
end

function prob=CDF(t,k,lambda)
%%CDF Evaluate the CDF of the Erlang probability distribution function.
%
%INPUTS: x The point(s) at which the Erlang CDF is to be evaluated.
%          Note that x>=0.
%        k The shape parameter of the Erlang distribution (k>=1).
%   lambda The rate parameter of the Erlang distribution (lambda>0).
%
%OUTPUTS: prob The value(s) of the Erlang CDF.
%
%The PDF of the Erlang distribution is given in Chapter 4 of [1] and also
%in [2]. The CDF can be obtained by integration.
%
%EXAMPLE:
%Here, we validate the CDF by generating random samples and comparing the
%empricial CDF with the CDF of this function.
% k=3;
% lambda=1/5;
% numSamples=1e5;
% numPoints=1000;
% x=linspace(0,60,numPoints);
% samples=ErlangD.rand([numSamples,1],k,lambda);
% CDFEmp=EmpiricalD.CDF(x,samples);
% CDF=ErlangD.CDF(x,k,lambda);
% 
% figure(1)
% clf
% plot(x,CDFEmp,'-k','linewidth',4)
% hold on
% plot(x,CDF,'-r','linewidth',2)
% h1=xlabel('x');
% h2=ylabel('CDF(x)');
% legend('Empirical CDF','CDF');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the empirical CDF matches well with the CDF of this
%function.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%[2] Weisstein, Eric W. "Erlang Distribution." From MathWorld--A Wolfram
%    Web Resource. https://mathworld.wolfram.com/ErlangDistribution.html
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    prob=gammainc(t*lambda,k);
end

function val=invCDF(prob,k,lambda)
%%INVCDF Evaluate the inverse of the cumulative distribution function of
%        the Erlang distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%        k The shape parameter of the Erlang distribution (k>=1).
%   lambda The rate parameter of the Erlang distribution (lambda>0).
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The CDF of the gamma distribution can be expressed in terms of the
%incomplete gamma function, as described in [1]. Since the Erlang
%distribution is a special case of the gamma distribution, the CDF of the
%Erlang distribution can be similarly expressed and the inverse incomplete
%gamma function built into Matlab can be used to invert it without the use
%of any toolboxes.
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% k=3;
% lambda=1/5;
% x=0.25;
% xBack=ErlangD.invCDF(ErlangD.CDF(x,k,lambda),k,lambda)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Gamma Distribution." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/GammaDistribution.html
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=gammaincinv(prob,k)/lambda;
end

function vals=rand(N,k,lambda)
%%RAND Generate Erlang distributed random variables with the specified
%      parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand 
%          returns an MXN1 matrix of random variables.
%        k The shape parameter of the Erlang distribution (k>=1).
%   lambda The rate parameter of the Erlang distribution (lambda>0).
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Erlang random variables.
%
%This is an implementation of the inverse transform algorithm of Chapter
%5.1 of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    U=rand(dims);
    vals=ErlangD.invCDF(U,k,lambda);
end

function val=momentGenFun(k,lambda,numDerivs,t)
%%MOMENTGENFUNCTION Evaluate the moment generating function (or one of its
%              derivatives) of the Erlang distribution. Taking the
%              numDerivth derivative of the moment generating function and
%              evaluating it at t=0 provides the numDerivth noncentral
%              moment of the distribution.
%
%INPUTS: k The shape parameter of the Erlang distribution (k>=1).
%   lambda The rate parameter of the Erlang distribution (lambda>0).
%  numDerivs The number of derivatives to take with respect to the
%          argument of the moment generating function. numDerivs>=0.
%        t The numPointsX1 or 1XnumPoints vector of points where the moment
%          generating function should be evaluated. If this parameter is
%          omitted or an empty matrix is passed, the default value of 0 is
%          used. This function only valid for t<lambda. 
%
%OUTPUTS: momentVal A numPointsX1 vector of the values of the derivatives
%                   of the moment generating function given at the points
%                   in t or at t=0 if t is omitted.
%
%The moment generating function of a random variable is defined to be
%E(exp(t*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter. It can be shown that the moment
%generating function of the Erlang distribution is
%f^0(t)=(1+t/lambda).^(-k);
%A general expression for the nth derivative can be found by identifying
%the pattern obtained from repeated differentiation.
%
%EXAMPLE:
%We verify the nth derivative of the meoment generating function evaluated
%at t=0 is the same as the output of nthMoment.
% n=12;
% k=3;
% lambda=1/4;
% theMoment=ErlangD.momentGenFun(k,lambda,n);
% theMomentAlt=ErlangD.nthMoment(n,k,lambda);
% RelErr=(theMoment-theMomentAlt)/theMoment
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(t))
       t=0; 
    end
    
    if(numDerivs==0)
        val=(1+t/lambda).^(-k);
    elseif(numDerivs==1)
        val=(k/lambda)*(1+t/lambda).^(-k-1);
    else
        val=abs(prod(-(1:(numDerivs-1))-k)*(k/lambda^numDerivs)*(1+t/lambda).^(-k-numDerivs));
    end
end

function val=entropy(k,lambda)
%%ENTROPY Obtain the differential entropy of the Erlang distribution given
%         in nats. The differential entropy of a continuous distribution is
%         entropy=-int_x p(x)*log(p(x)) dx where the integral is over all
%         values of x. Units of nats mean that the natural logarithm is
%         used in the definition. Unlike the Shannon entropy for discrete
%         variables, the differential entropy of continuous variables can
%         be both positive and negative.
%
%INPUTS: k The shape parameter of the Erlang distribution (k>=1).
%   lambda The rate parameter of the Erlang distribution (lambda>0).
%
%OUTPUTS: entropyVal The value of the differential entropy in nats.
%
%Differential entropy is defined in Chapter 8 of [1].
%
%REFERENCES:
%[1] T. M. Cover and J. A. Thomas, Elements of Information Theory, 2nd ed.
%    Hoboken, NJ: Wiley-Interscience, 2006.
%
%September 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=k+log(gamma(k)/lambda)-(k-1)*polygamma(0,k);
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
