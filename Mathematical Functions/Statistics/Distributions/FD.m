classdef FD
%%FD Functions to handle the scalar F-distribution. If x is a random
%    variable with a central chi-squared distribution with m degrees of
%    freedom and y is a central chi squared random variable with n degrees
%    of freedom, then as noted in Chapter 6-3 of [1], the ratio (x/m)/(y/n)
%    is given by the F distribution with (m,n) degrees of freedom. This is
%    also known as the Fisher-Snedecor distribution and as Fisher's
%    z-distribution.
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand, kthMoment
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=mean(n)
%%MEAN Obtain the mean of F probability distribution.
% 
%INPUTS: n The number of degrees of freedom of the denominator term forming
%          the F distribution. n>2 for the mean to exist.
%
%OUTPUTS: val The mean of the F distribution under consideration.
%
%The mean on the page opposite the inner front cover of [1] and in [2]. The
%mean only exists for n>2.
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% m=3;
% n=7;
% meanVal=FD.mean(n)
% numSamp=1e6;
% sampMeanVal=mean(FD.rand([numSamp,1],m,n))
%One will see that both values are about 1.4.
% 
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%[2] Weisstein, Eric W. "Fisher's z-Distribution." From MathWorld--A
%    Wolfram Web Resource.
%    http://mathworld.wolfram.com/Fishersz-Distribution.html
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(n<=2)
        error('The mean does not exist for n<=2');
    end

    val=n/(n-2);
end
    
function val=var(m,n)
%%VAR Obtain the variance of F probability distribution.
%
%INPUTS: m The number of degrees of freedom of the numerator term forming
%          the F distribution.
%        n The number of degrees of freedom of the denominator term forming
%          the F distribution.
%
%OUTPUTS: val The variance of the F distribution under consideration.
%
%The variance is given on the page opposite the inner front cover of [1].
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% m=3;
% n=7;
% varVal=FD.var(m,n)
% numSamp=1e6;
% sampVarVal=var(FD.rand([numSamp,1],m,n))
%One will get both values around 3.48, but there can be a large amount of
%variation between runs than one might expect.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(n<=4)
        error('The variance does not exist for n<=4');
    end

    val=(2*n^2*(m+n-2))/(m*(n-4)*(n-2)^2);
end
    
function val=PDF(x,m,n)
%%PDF Evaluate the F probability distribution function (PDF) at the desired
%     points.
%
%INPUTS: x The point(s) at which the F PDF is to be evaluated. Note that
%          x>=0 for nonzero probabilities.
%        m The number of degrees of freedom of the numerator term forming
%          the F distribution.
%        n The number of degrees of freedom of the denominator term forming
%          the F distribution.
%
%OUTPUTS: val The value(s) of the F PDF evaluated at x.
%
%The PDF is given in Chapter 6-3 of [1] and in [2].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% m=3;
% n=7;
% numSamples=1e5;
% figure(1)
% clf
% histogram(FD.rand([numSamples,1],m,n),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,20,numPoints);
% vals=FD.PDF(x,m,n);
% plot(x,vals,'linewidth',2)
% axis([0,10,0,0.7])
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
%[2] Weisstein, Eric W. "Fisher's z-Distribution." From MathWorld--A
%    Wolfram Web Resource.
%    http://mathworld.wolfram.com/Fishersz-Distribution.html
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=(m/n)*exp((m/2-1)*log(m/n*x)-((m+n)/2)*log(1+(m/n)*x)-betaln(m/2,n/2));
    val(x<=0)=0;
end
    
function prob=CDF(x,m,n)
%%CDF Evaluate the cumulative distribution function (CDF) of the F
%     distribution at desired points.
%
%INPUTS: x The point(s) at which the F CDF is to be evaluated.
%          Note that x>0 for nonzero values.
%        m The number of degrees of freedom of the numerator term forming
%          the F distribution.
%        n The number of degrees of freedom of the denominator term forming
%          the F distribution.
%
%OUTPUTS: prob The value(s) of the CDF of the F distribution evaluated at
%              x.
%
%The PDF is given in Chapter 6-3 of [1] and in [2]. Integrating the PDF can
%be shown to give an incomplete regularized beta function.
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% m=3;
% n=7;
% x=3;
% numSamples=1e5;
% prob=FD.CDF(x,m,n)
% probSamp=mean(FD.rand([numSamples,1],m,n)<=x)
%One will see that both values are about 0.895.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%[2] Weisstein, Eric W. "Fisher's z-Distribution." From MathWorld--A
%    Wolfram Web Resource.
%    http://mathworld.wolfram.com/Fishersz-Distribution.html
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    prob=betainc(m*x./(m*x+n),m/2,n/2);
end

function val=invCDF(prob,m,n)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the F distribution. 
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%           m The number of degrees of freedom of the numerator term
%             forming the F distribution.
%           n The number of degrees of freedom of the denominator term
%             forming the F distribution.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The PDF is given in Chapter 6-3 of [1] and in [2]. Integrating the PDF can
%be shown to give an incomplete regularized beta function. The betaincinv
%function is used to take the inverse.
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% m=3;
% n=7;
% x=3;
% xBack=FD.invCDF(FD.CDF(x,m,n),m,n)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%[2] Weisstein, Eric W. "Fisher's z-Distribution." From MathWorld--A
%    Wolfram Web Resource.
%    http://mathworld.wolfram.com/Fishersz-Distribution.html
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    temp=betaincinv(prob,m/2,n/2);
    val=(n*temp)./(m*(1-temp));
end

    
function vals=rand(N,m,n)
%%RAND Generate F distributed random variables with the given parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%        m The number of degrees of freedom of the numerator term forming
%          the F distribution.
%        n The number of degrees of freedom of the denominator term forming
%          the F distribution.
%
%%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated F random variables.
%
%The random variables are obtained by dividing two central chi squared
%distributions as discussed in Chapter 6-3 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    vals=(ChiSquareD.rand(dims,m)/m)./(ChiSquareD.rand(dims,n)/n);
end

function val=kthMoment(k,m,n)
%%KTHMOMENT Evaluate the kth noncentral moment of the F distribution, if it
%           exist, or throw an error if it does not exist.
%
%INPUTS: k The positive integer order of the moment.
%        m The number of degrees of freedom of the numerator term forming
%          the F distribution.
%        n The number of degrees of freedom of the denominator term forming
%          the F distribution.
%
%OUTPUTS: val The value of the kth moment.
%
%The PDF is given in Chapter 6-3 of [1] and in [2].The kth moment can be
%evaluated symbolically, for example in Mathematica. The resulting moment
%only exists for n>2*k.
%
%EXAMPLE:
%Here, we validate the kthMoment by comparing it to a sample moment.
% m=3;
% n=7;
% k=2;
% kMoment=FD.kthMoment(k,m,n)
% numSamples=1e6;
% kMomentSample=mean(FD.rand([numSamples,1],m,n).^k)
%One will find both are about 5.4. Note that the variance of the sample
%moment increases greatly as k increases.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%[2] Weisstein, Eric W. "Fisher's z-Distribution." From MathWorld--A
%    Wolfram Web Resource.
%    http://mathworld.wolfram.com/Fishersz-Distribution.html
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(n<=2*k)
       error('The moment does not exist for the given parameters.') 
    end
    
    val=exp(k.*log(n./m)+gammaln(k+m/2)+gammaln(n/2-k)-gammaln(m/2)-gammaln(n/2));
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
