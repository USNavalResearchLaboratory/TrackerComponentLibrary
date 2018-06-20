classdef GumbelD
%%GUMBELD Functions to handle the Gumbel distribution. This is also known
%         as the generalized extreme value distribution type-I, and the
%         log-Weibull distribution.
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.    
    
methods(Static)
function val=mean(mu,betaVal)
%%MEAN Obtain the mean of the Gumbel distribution.
%
%INPUTS: mu The real location parameter of the distribution 
%   betaVal The scale parameter of the distribution. betaVal>0.
%
%OUTPUTS: val The mean of the Gumbel distribution.
%
%The Gumbel distribution is described in [1].
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% mu=-5;
% betaVal=3;
% meanVal=GumbelD.mean(mu,betaVal)
% numSamp=1e6;
% sampMeanVal=mean(GumbelD.rand([numSamp,1],mu,betaVal))
%One will find both values are about -3.268
%
%REFERENCES:
%[1] Weisstein, Eric W. "Gumbel Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/GumbelDistribution.html
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    gammaVal=-psi(1);%The Euler-Mascheroni constant.
    val=mu+betaVal*gammaVal;
end

function val=var(betaVal)
%%VAR Obtain the variance of the Gumbel distribution.
%
%INPUTS: betaVal The scale parameter of the distribution. betaVal>0.
%
%OUTPUTS: val The variance of the Gumbel distribution.
%
%The Gumbel distribution is described in [1].
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% mu=-5;
% betaVal=3;
% varVal=GumbelD.var(betaVal)
% numSamp=1e6;
% sampVarVal=var(GumbelD.rand([numSamp,1],mu,betaVal))
%One will find both values are about 14.804
%
%REFERENCES:
%[1] Weisstein, Eric W. "Gumbel Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/GumbelDistribution.html
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

   val=pi^2/6*betaVal^2;
end

function val=PDF(x,mu,betaVal)
%%PDF Evaluate the probability density function (PDF) of the Gumbel
%     distribution at one or more points.
%
%INPUTS: x The point(s) at which the Gumbel PDF is to be evaluated.
%       mu The real location parameter of the distribution 
%  betaVal The scale parameter of the distribution. betaVal>0.
%
%OUTPUTS: val The value(s) of the Gumbel PDF.
%
%The Gumbel distribution is described in [1].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% mu=-5;
% betaVal=3;
% numSamples=1e5;
% figure(1)
% clf
% histogram(GumbelD.rand([numSamples,1],mu,betaVal),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(-15,15,numPoints);
% vals=GumbelD.PDF(x,mu,betaVal);
% plot(x,vals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Gumbel Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/GumbelDistribution.html
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    z=(x-mu)/betaVal;
    val=(1/betaVal)*exp(-(z+exp(-z)));
end

function prob=CDF(x,mu,betaVal)
%%CDF Evaluate the cumulative distribution function (CDF) of the gumbel
%     distribution at one or more desired points.
%
%INPUTS: x The point(s) at which the Gumbel CDF is to be evaluated.
%       mu The real location parameter of the distribution 
%  betaVal The scale parameter of the distribution. betaVal>0.
%
%OUTPUTS: prob The value(s) of the Gumbel CDF.
%
%The Gumbel distribution is described in [1].
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% mu=-5;
% betaVal=3;
% x=-1;
% numSamples=1e5;
% prob=GumbelD.CDF(x,mu,betaVal)
% probSamp=mean(GumbelD.rand([numSamples,1],mu,betaVal)<=x)
%One will see that both values are about 0.768.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Gumbel Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/GumbelDistribution.html
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    prob=exp(-exp(-(x-mu)/betaVal));
end
    
function x=invCDF(prob,mu,betaVal)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the Gumbel distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%       mu The real location parameter of the distribution 
%  betaVal The scale parameter of the distribution. betaVal>0.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The CDF of the Gumbel distribution can be easily algebraicly inverted,
%which is what is done here. The Gumbel distribution is described in [1].
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% mu=-5;
% betaVal=3;
% x=-1;
% xBack=GumbelD.invCDF(GumbelD.CDF(x,mu,betaVal),mu,betaVal)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Gumbel Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/GumbelDistribution.html
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    x=mu+betaVal*log(1./log(1./prob));
end

function x=rand(N,mu,betaVal)
%%RAND Generate Gumbel distributed random variables.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand 
%          returns an MXN1 matrix of random variables.
%       mu The real location parameter of the distribution 
%   betaVal The scale parameter of the distribution. betaVal>0.
%
%OUTPUTS: x A matrix whose dimensions are determined by N of the generated
%           Gumbel random variables.
%
%This is an implementation of the inverse transform algorithm of Chapter
%5.1 of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    U=rand(dims);

    x=GumbelD.invCDF(U,mu,betaVal);
end

function entropyVal=entropy(betaVal)
%%ENTROPY Obtain the differential entropy of the Gumbel distribution given
%         in nats. The differential entropy of a continuous distribution is
%         entropy=-int_x p(x)*log(p(x)) dx where the integral is over all
%         values of x. Units of nats mean that the natural logarithm is
%         used in the definition. Unlike the Shannon entropy for discrete
%         variables, the differential entropy of continuous variables can
%         be both positive and negative.
%
%INPUTS: betaVal The scale parameter of the distribution. betaVal>0.
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

    gammaVal=-psi(1);%The Euler-Mascheroni constant.
    entropyVal=log(betaVal)+gammaVal+1;
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
