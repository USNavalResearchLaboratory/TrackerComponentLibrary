classdef RayleighD
%%RAYLEIGHD Functions to handle the Rayleigh distribution.
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand, momentGenFun,
%                         entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=mean(sigma)
%%MEAN Obtain the mean of the Rayleigh distribution for a given parameter.
%
%INPUTS: sigma The parameter of the Rayleigh distribution.
%
%OUTPUTS: val The mean of the Rayleigh distribution.
%
%The mean of the Rayleigh distribution is given on the inside cover of [1].
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% sigma=4;
% meanVal=RayleighD.mean(sigma)
% numSamp=1e6;
% sampMeanVal=mean(RayleighD.rand([numSamp,1],sigma))
%One will find both values are about 5.013
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    val=sigma*sqrt(pi/2);
end

function val=var(sigma)
%%VAR Obtain the variance of the Rayleigh distribution for a given
%     parameter.
%
%INPUTS: sigma The parameter of the Rayleigh distribution.
%
%OUTPUTS:  val The variance of the Rayleigh distribution.
%
%The variance of the Rayligh distribution is given on the inside cover of
%[1].
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% sigma=4;
% varVal=RayleighD.var(sigma)
% numSamp=1e6;
% sampVarVal=var(RayleighD.rand([numSamp,1],sigma))
%One will find both values are about 6.867
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    val=(2-pi/2)*sigma^2;
end  
    
function val=PDF(x,sigma)
%%PDF Evaluate the Rayleigh probability distribution function (PDF) at one
%     or more desired points.
%
%INPUTS: x The point(s) at which the Rayleigh PDF is to be  evaluated.
%    sigma The parameter of the Rayleigh distribution.
%
%OUTPUTS: val The value(s) of the Rayleigh PDF with parameter sigma
%             evaluated at x.
%
%The PDF of the Rayleigh distribution is given in Chapter 4.3 of [1].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% sigma=5;
% numSamples=1e5;
% figure(1)
% clf
% histogram(RayleighD.rand([numSamples,1],sigma),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,20,numPoints);
% vals=RayleighD.PDF(x,sigma);
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
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(x));
    sel=(x>=0);

    val(sel)=(x(sel)/sigma^2).*exp(-x(sel).^2/(2*sigma^2));
end

function val=CDF(x,sigma)
%%CDF Evaluate the cumulative distribution function (CDF) of the Rayleigh
%     distribution at one or more desired points.
%
%INPUTS: x The point(s) at which the Rayleigh CDF is to be evaluated.
%    sigma The parameter of the Rayleigh distribution.
%
%OUTPUTS: val The value(s) of the Rayleigh CDF with parameter sigma
%             evaluated at x.
%
%The CDF of the Rayleigh distribution can be easily derived by integration
%the PDF of the Rayleigh distribution, which is given in Chapter 4.3 of
%[1].
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% sigma=5;
% x=3;
% numSamples=1e5;
% prob=RayleighD.CDF(x,sigma)
% probSamp=mean(RayleighD.rand([numSamples,1],sigma)<=x)
%One will see that both values are about 0.1647.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(x));
    sel=(x>=0);

    val(sel)=1-exp(-x(sel).^2/(2*sigma^2));
end

function val=invCDF(prob,sigma)
%%INVCDF Evaluate the inverse of the cumulative distribution function of
%        the Rayleigh distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%       sigma The parameter of the Rayleigh distribution.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The CDF of the Rayleigh distribution is simple and can be easily inverted
%using logarithms. The CDF of the Rayleigh distribution can be easily
%derived by integration the PDF of the Rayleigh distribution, which is
%given in Chapter 4.3 of [1].
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% sigma=5;
% x=3;
% xBack=RayleighD.invCDF(RayleighD.CDF(x,sigma),sigma)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=sigma*sqrt(-2*log(1-prob));
end

function vals=rand(N,sigma)
%%RAND Generate Rayleigh random variables with a given parameter sigma.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of  random variables.
%    sigma The parameter of the Rayleigh distribution.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Rayleigh random variables.
%
%This generates Rayleigh distributed random variables by transforming 
%normally distributed random variables using the identity given in
%Chapter 2.1.4 of [1].
%
%REFERENCES:
%[1] J. G. Proakis, Digital Communication. Ed. 4, Boston, MA: McGraw Hill,
%    2001.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N N];
    else
        dims=N;
    end

    X=sigma*randn(dims);
    Y=sigma*randn(dims);

    vals=sqrt(X.^2+Y.^2);
end

function momentVal=momentGenFun(sigma,numDerivs,t)
%%MOMENTGENFUN Evaluate the moment generating function (or one of its
%              derivatives) of the exponential distribution. Taking the kth
%              derivative of the moment generating function and evaluating
%              it at t=0 provides the kth noncentral moment of the
%              distribution.
%
%INPUTS:  sigma The parameter of the Rayleigh distribution. 
%     numDerivs The number of derivatives to take with respect to the
%               argument of the moment generating function. numDerivs>=0.
%             t The numPointsX1 or 1XnumPoints vector of points where the
%               moment generating function should be evaluated. If this
%               parameter is omitted or an empty matrix is passed, the
%               default value of 0 is used.
%
%OUTPUTS: momentVal A numPointsX1 vector of the values of the derivatives
%                   of the moment generating function given at the points
%                   in t or at t=0 if t is omitted.
%
%The moment generating function of a random variable is defined to be
%E(exp(t*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter. It can be shown that the moment
%generating function of the Rayleigh distribution is
%f^0(t)=1+sigma*t*exp(sigma^2*t^2/2)*sqrt(pi/2)*(erf(sigma*t/sqrt(2))+1);
%This function finds derivatives of the moment generating function by
%realizing that the derivatives of the moment generating function consists
%of three parts: A polynomial part, a polynomial times
%exp(sigma^2*t^2/2)*erf(sigma*t/sqrt(2)) and a polynomial times
%exp(sigma^2*t^2/2). Differentiation is thus performed by just keeping
%track of those three polynomials.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(nargin<3||isempty(t))
       t=0; 
    end
    
    %If here, we must take at least one derivative. We shall keep track of
    %the polynomials times each type of coefficient. Polynomial
    %coefficients are ordered in the order that polyder and polyval take:
    %highest order coefficient first. 
    %poly0 is the polynomial term that is not multiplied by anything.
    poly0=1;
    %polyExpErf is the polynomial term that is multiplied by
    %exp(sigma^2*t^2/2)*erf(sigma*t/sqrt(2))
    polyExpErf=[sigma*sqrt(pi/2),0];
    %polyExp is the polynomial term that is multiplied by
    %exp(sigma^2*t^2/2)
    polyExp=[sigma*sqrt(pi/2),0];
    
    %Take derivatives.
    while(numDerivs>0)
        %%%%Term 1:
        %The derivative of the polynomial-only portion is straightforward.
        poly0=polyder(poly0);
        
        %%%%Term 2:
        %For the derivative of the erf-exponent product, we take two steps
        %due to the chain rule. First, differentiate the polynomial part.
        %That will remain multiplied by the erf-exponent product.
        polyExpErfDeriv=polyder(polyExpErf);
        %Next, we differentiate the non-polynomial part. This results in a
        %constant term that is multiplied by the polyExpErf and must be
        %added to the polynomial-only portion, as well as a linear term
        %affecting the erf-exponent product.
        
        %Deal with the polynomial term first.
        term2Add=sigma*sqrt(2/pi)*polyExpErf;
        term2AddLen=length(term2Add);
        poly0Len=length(poly0);
        if(term2AddLen<poly0Len)
            poly0((end-term2AddLen+1):end)=poly0((end-term2AddLen+1):end)+term2Add;
        else
            term2Add((end-poly0Len+1):end)=term2Add((end-poly0Len+1):end)+poly0;
            poly0=term2Add;
        end
        
        %Next, deal with the part times the erf-exponent product.
        term2Add=conv([sigma^2,0],polyExpErf);%Polynomial multiplication
        term2AddLen=length(term2Add);
        polyExpErfDerivLen=length(polyExpErfDeriv);
        if(term2AddLen<polyExpErfDerivLen)
            polyExpErfDeriv((end-term2AddLen+1):end)=polyExpErfDeriv((end-term2AddLen+1):end)+term2Add;
        else
            term2Add((end-polyExpErfDerivLen+1):end)=term2Add((end-polyExpErfDerivLen+1):end)+polyExpErfDeriv;
            polyExpErfDeriv=term2Add;
        end
        polyExpErf=polyExpErfDeriv;
        
        %%Term 3:
        %The derivative of the exp product is similar to that of the
        %exp-erf product, except there is no new polynomial-only term.
        %The chain rule is used again. First, we differentiate the
        %polynomial part.
        polyExpDeriv=polyder(polyExp);
        %Next, we differentiate the non-polynomial part and add it to the
        %polynomial part.
        term2Add=conv([sigma^2,0],polyExp);%polynomial multiplication
        term2AddLen=length(term2Add);
        polyExpDerivLen=length(polyExpDeriv);
        if(term2AddLen<polyExpErfDerivLen)
            polyExpDeriv((end-term2AddLen+1):end)=polyExpDeriv((end-term2AddLen+1):end)+term2Add;
        else
            term2Add((end-polyExpDerivLen+1):end)=term2Add((end-polyExpDerivLen+1):end)+polyExpDeriv;
            polyExpDeriv=term2Add;
        end
        polyExp=polyExpDeriv;
        numDerivs=numDerivs-1;
    end
    
    %If here, all of the derivatives have been taken and the function can
    %be evaluated.
    momentVal=polyval(poly0,t)+polyval(polyExpErf,t)*exp(sigma^2*t^2/2)*erf(sigma*t/sqrt(2))+polyval(polyExp,t)*exp(sigma^2*t^2/2);
end

function entropyVal=entropy(sigma)
%%ENTROPY Obtain the differential entropy of the Rayleigh distribution
%         given in nats. The differential entropy of a continuous
%         distribution is entropy=-int_x p(x)*log(p(x)) dx where the
%         integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: sigma The parameter of the Rayleigh distribution.
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
    entropyVal=1+log(sigma/sqrt(2))+gammaVal/2;
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
