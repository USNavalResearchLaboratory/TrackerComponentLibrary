classdef WeibullD
%%WEIBULLD Functions to handle the Weibull distribution.
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=mean(lambda,k)
%%MEAN Obtain the mean of the Weibull distribution for given scale and
%      shape parameters.
%
%INPUTS: lambda The scale parameter of the distribution. lambda>0.
%             k The shape parameter of the distribution. k>0.
%
%OUTPUTS:val The mean of the Weibull distribution.
%
%The mean of the Weibull distribution is given on the inside of the front
%cover of [1].
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% lambda=2.2;
% k=3.1;
% meanVal=WeibullD.mean(lambda,k)
% numSamp=1e6;
% sampMeanVal=mean(WeibullD.rand([numSamp,1],lambda,k))
%One will find both values are about 1.967
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

   val=lambda*gamma(1+1/k); 
end

function val=var(lambda,k)
%%VAR Obtain the variance of the Weibull distribution for given scale and
%     shape parameters/
%
%INPUTS: lambda The scale parameter of the distribution. lambda>0.
%             k The shape parameter of the distribution. k>0.
%
%OUTPUTS:val The variance of the Weibull distribution.   
%
%The variance of the Weibull distribution is given on the inside of the front
%cover of [1].
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% lambda=2.2;
% k=3.1;
% varVal=WeibullD.var(lambda,k)
% numSamp=1e6;
% sampVarVal=var(WeibullD.rand([numSamp,1],lambda,k))
%One will find both values are about 0.4821
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%    
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

   val=lambda^2*(gamma(1+2/k)-(gamma(1+1/k))^2); 
end

function val=PDF(x,lambda,k)
%%PDF Evaluate the probability density function (PDF) of the Weibull
%     distribution at one or more points.
%
%INPUTS: x The point(s) at which the Weibull PDF is to be evaluated. x>=0
%          for nonzero PDF values.
%   lambda The scale parameter of the distribution. lambda>0.
%        k The shape parameter of the distribution. k>0.
%
%OUTPUTS: val The value(s) of the Weibull PDF.
%
%The PDF of the Weibull distribution is given on the inside of the front
%cover of [1].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% lambda=2.2;
% k=3.1;
% numSamples=1e5;
% figure(1)
% clf
% histogram(WeibullD.rand([numSamples,1],lambda,k),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,6,numPoints);
% vals=WeibullD.PDF(x,lambda,k);
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
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(k/lambda)*(x/lambda).^(k-1).*exp(-(x/lambda).^k);
    val(x<0)=0;
end

function val=CDF(x,lambda,k)
%%CDF Evaluate the cumulative distribution function (CDF) of the Weibull
%     distribution at one or more desired points.
%
%INPUTS: x The point(s) at which the Weibull CDF is to be evaluated.
%   lambda The scale parameter of the distribution. lambda>0.
%        k The shape parameter of the distribution. k>0.
%
%OUTPUTS: val The value(s) of the Weibull CDF.
%
%The CDF is the integral from 0 to x of the PDF. The Weibull PDF is
%fairely simple and the integral is not hard. The Weibull PDF is given on
%the inside of the front cover of [1].
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% lambda=2.2;
% k=3.1;
% x=2;
% numSamples=1e5;
% prob=WeibullD.CDF(x,lambda,k)
% probSamp=mean(WeibullD.rand([numSamples,1],lambda,k)<=x)
%One will see that both values are about 0.5249.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=1-exp(-(x/lambda).^k);
    val(x<0)=0;
end

function val=invCDF(prob,lambda,k)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the Weibull distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%      lambda The scale parameter of the distribution. lambda>0.
%           k The shape parameter of the distribution. k>0.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The CDF of the Weibull distribution can be easily algebraicly inverted,
%which is what is done here. The Weibull CDF is the integral from 0 to x of
%the Weibull PDF, which is given on the inside of the front cover of [1].
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% lambda=2.2;
% k=3.1;
% x=2;
% xBack=WeibullD.invCDF(WeibullD.CDF(x,lambda,k),lambda,k)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=lambda*(-log(1-prob)).^(1/k);

end
    
function vals=rand(N,lambda,k,lowerBound)
%%RAND Generate Weibull distributed random variables with the given
%      parameters. This function can also be used to generate random
%      variable from a Weibull distribution that is clipped so that only
%      values above a certain threshold are generated. That can be useful
%      when performing simulations and generating false alarms that are all
%      above a given threshold.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%   lambda The scale parameter of the distribution. lambda>0.
%        k The shape parameter of the distribution. k>0.
% lowerBound An optional parameter specifying a lower bound below which
%          none of the randomly generated values should go. If this is
%          omitted or an empty matrix is passed, then lowerBound=0 and a
%          non-clipped Weibull distribution is used.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Weibull random variables.
%
%This is an implementation of the inverse transform algorithm of Chapter
%5.1 of [1]. To clip values to only being above lowerBound, the CDF value
%of the lower bound is determined and the uniform random variables as used
%in the inverse transform algorithm are all only generated above the
%threshold.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%  
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    if(nargin<4||isempty(lowerBound)||lowerBound<=0)
        %If we are evaluating a standard Weibull-distributed random
        %variable rather than one that has been clipped to a particular
        %region.
        U=rand(dims);
        vals=WeibullD.invCDF(U,lambda,k);
    else
        %Only generate Weibull random variables that are >=lowerBound.
        p=WeibullD.CDF(lowerBound,lambda,k);
        %Generate a uniform random variable between p and 1.
        U=p+(1-p)*rand(dims);
        vals=WeibullD.invCDF(U,lambda,k);
    end
end

function entropyVal=entropy(lambda,k)
%%ENTROPY Obtain the differential entropy of the Weibull distribution
%         given in nats. The differential entropy of a continuous
%         distribution is entropy=-int_x p(x)*log(p(x)) dx where the
%         integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: lambda The scale parameter of the distribution. lambda>0.
%        k The shape parameter of the distribution. k>0.
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

    gammaVal=-psi(1);
    entropyVal=gammaVal*(1-1/k)+log(lambda/k)+1;
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
