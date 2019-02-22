classdef CauchyD
%%CAUCHYD Functions to handle the scalar Cauchy distribution. For
%         multivariate distributions, use StudentTD with one degree of
%         freedom. The Cauchy distribution has very long tails and no
%         defined mean. In [1], it was noted to arise in various statistics
%         related to directional measurements in the presence of glint. The
%         Cauchy distribution is stable, which means that the sum of two
%         Cauchy random variables is also a Cauchy random variable.
%Implemented methods are: PDF, CDF, invCDF, rand, entropy
%
%REFERENCES:
%[1] U. Nickel, "Angular superresolution with phased array radar: A review
%    of algorithms and operational constraints," IEEE Proceedings, vol.
%    134, Part F, no. 1, pp. 53-59, Feb. 1987.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=PDF(x,mu,alphaVal)
%%PDF Evaluate the probability density function (PDF) of the scalar Cauchy
%     distribution.
%
%INPUTS: x The point or a matrix of points at which the PDF of the Cauchy
%          distribution should be evaluated.
%       mu The location parameter of the distribution, which is also its
%          mode.
% alphaVal The scale parameter of the distribution.
%
%OUTPUTS: val The value(s) of the PDF evaluated at the point(s) in x.
%
%The PDF is given in Chapter 4 of [1].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% mu=4;
% alphaVal=2;
% numSamples=1e6;
% figure(1)
% clf
% histogram(CauchyD.rand([numSamples,1],mu,alphaVal),'Normalization','pdf','BinLimits',[-50,50])
% hold on
% numPoints=1000;
% x=linspace(-20,20,numPoints);
% vals=CauchyD.PDF(x,mu,alphaVal);
% plot(x,vals,'linewidth',2)
% axis([-20,20,0,0.2])
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
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    val=1./(pi*alphaVal*(1+((x-mu)/alphaVal).^2));
end

function prob=CDF(x,mu,alphaVal)
%%CDF Evaluate the cumulative distribution function (CDF) of the scalar
%     Cauchy distribution.
%
%INPUTS: x The point or a matrix of points at which the CDF of the Cauchy
%          distribution should be evaluated.
%       mu The location parameter of the distribution, which is also its
%          mode.
% alphaVal The scale parameter of the distribution.
%
%OUTPUTS: prob The value(s) of the CDF of the Cauchy distribution evaluated
%              at x.
%
%The CDF of the Cauchy distribution is given in [1].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% mu=4;
% alphaVal=2;
% x=5;
% numSamples=1e6;
% prob=CauchyD.CDF(x,mu,alphaVal)
% probSamp=mean(CauchyD.rand([numSamples,1],mu,alphaVal)<=x)
%One will see that both values are about 0.6476.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Cauchy Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/CauchyDistribution.html
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    prob=(1/pi)*atan((x-mu)/alphaVal)+1/2;
end

function x=invCDF(prob,mu,alphaVal)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the scalar Cauchy distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%          mu The location parameter of the distribution, which is also its
%             mode.
%    alphaVal The scale parameter of the distribution.
%
%OUTPUTS: x The argument(s) of the CDF that would give the probability or
%           probabilities in prob.
%
%The CDF of the Cauchy distribution is given in [1]. its inversion is
%straightforward.
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% mu=4;
% alphaVal=2;
% x=5;
% xBack=CauchyD.invCDF(CauchyD.CDF(x,mu,alphaVal),mu,alphaVal)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Cauchy Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/CauchyDistribution.html
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    x=alphaVal*tan(pi*(prob-1/2))+mu;
end

function x=rand(N,mu,alphaVal)
%%RAND Generate scalar Cauchy distributed random variables.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%       mu The location parameter of the distribution, which is also its
%          mode.
% alphaVal The scale parameter of the distribution.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Cauchy random variables.
%
%This function implements the inverse transformation method of Chapter 5.1
%of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    if(isscalar(N))
        dims=[N N];
    else
        dims=N;
    end
    
    U=rand(dims);
    x=CauchyD.invCDF(U,mu,alphaVal);
end

function entropyVal=entropy(alphaVal)
%%ENTROPY Obtain the differential entropy of the Cauchy distribution given
%         in nats. The differential entropy of a continuous distribution is
%         entropy=-int_x p(x)*log(p(x))  dx where the integral is over all
%         values of x. Units of nats mean that the natural logarithm is
%         used in the definition. Unlike the Shannon entropy for discrete
%         variables, the differential entropy of continuous variables can
%         be both positive and negative.
%
%INPUTS: alphaVal The scale parameter of the distribution.
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
    
    entropyVal=log(4*pi*alphaVal);
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
