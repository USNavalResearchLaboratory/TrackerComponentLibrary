classdef ExponentialD
%%EXPONENTIALD Functions to handle the exponential distribution.
%Implemented methods are: mean, var, PDF, CDF, invCDF, momentGenFun,
%                         cumGenFun, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static) 
function val=mean(lambda)
%%MEAN Obtain the mean of the exponential distribution for a given rate
%      parameter.
%
%INPUTS: lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS: val The mean of the exponential distribution.
%
%The mean of the exponential distribution is given in Chapter 2.9 of [1].
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% lambda=4;
% meanVal=ExponentialD.mean(lambda)
% numSamp=1e6;
% sampMeanVal=mean(ExponentialD.rand([numSamp,1],lambda))
%One will find both values are about 0.25
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.    
    
    val=1/lambda;
end

function val=var(lambda)
%%VAR Obtain the variance of the exponential distribution for a given rate
%     parameter.
%
%INPUTS: lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS: val The variance of the exponential distribution.
%
%The variance of the exponential distribution is given in Chapter 2.9 of
%[1].
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% lambda=4;
% varVal=ExponentialD.var(lambda)
% numSamp=1e6;
% sampVarVal=var(ExponentialD.rand([numSamp,1],lambda))
%One will find both values are about 0.0625
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.    
    
    val=1/lambda^2;
end

function val=PDF(x,lambda)
%%PDF Evaluate the PDF of the exponential probability distribution function
%     at one or more desired points.
%
%INPUTS: x The point(s) at which the expoenntial distribution is to be 
%          evaluated.
%   lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS: val The value(s) of the exponential PDF with given rate parameter
%             evaluated at x.
%
%The PDF of the exponential distribution is given in Chapter 2.9 of [1].
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% lambda=5;
% numSamples=1e5;
% figure(1)
% clf
% histogram(ExponentialD.rand([numSamples,1],lambda),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,2,numPoints);
% vals=ExponentialD.PDF(x,lambda);
% plot(x,vals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(x));
    sel=(x>=0);
    
    val(sel)=lambda*exp(-lambda*x(sel));
end

function val=CDF(x,lambda)
%%CDF Evaluate the CDF of the exponential probability distribution function
%     at one or more desired points.
%
%INPUTS: x The point(s) at which the exponential CDF is to be evaluated.
%          Note that x>=0.
%   lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS: prob The value(s) of the expoenntial CDF with given rate
%              parameter evaluated at x.
%
%The CDF of the exponential distribution is given in Chapter 2.9 of [1].
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% lambda=5;
% x=0.25;
% numSamples=1e5;
% prob=ExponentialD.CDF(x,lambda)
% probSamp=mean(ExponentialD.rand([numSamples,1],lambda)<=x)
%One will see that both values are about 0.713.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.   
    
    val=zeros(size(x));
    sel=(x>=0);
    
    val(sel)=1-exp(-lambda*x(sel));
end

function val=invCDF(prob,lambda)
%%INVCDF Evaluate the inverse of the cumulative distribution function of
%        the exponential distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%      lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The CDF of the exponential distribution is very simple and is easily
%inverted using logarithms, as is done here. The CDF of the exponential
%distribution is given in Chapter 2.9 of [1].
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% lambda=5;
% x=0.25;
% xBack=ExponentialD.invCDF(ExponentialD.CDF(x,lambda),lambda)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.    
    
    val=-log(1-prob)/lambda;
end

function momentVal=momentGenFun(lambda,numDerivs,t)
%%MOMENTGENFUN Evaluate the moment generating function (or one of its
%              derivatives) of the exponential distribution. Taking the kth
%              derivative of the moment generating function and evaluating
%              it at t=0 provides the kth noncentral moment of the
%              distribution.
%
%INPUTS: lambda The rate parameter of the distribution. lambda>0.
%     numDerivs The number of derivatives to take with respect to the
%               argument of the moment generating function. numDerivs>=0.
%             t The numPointsX1 or 1XnumPoints vector of points where the
%               moment generating function should be evaluated. If this
%               parameter is omitted or an empty matrix is passed, the
%               default value of 0 is used. This function only valid for
%               t<lambda.         
%
%OUTPUTS: momentVal A numPointsX1 vector of the values of the derivatives
%                   of the moment generating function given at the points
%                   in t or at t=0 if t is omitted.
%
%The moment generating function of a random variable is defined to be
%E(exp(t*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter. It can be shown that the moment
%generating function of the exponential distribution is
%f^0(t)=lambda/(lambda-t)
%A general expression for the nth derivative can be found by identifying
%the pattern obtained from repeated differentiation.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(t))
       t=0; 
    end
    
    momentVal=factorial(numDerivs)*lambda/(lambda-t)^(numDerivs+1);
end

function cumVal=cumGenFun(lambda,numDerivs,t)
%%CUMGENFUN Evaluate the cumulant generating function (or one of its
%           derivatives) of the exponential distribution. The cumulant
%           generating function is the natural logarithm of the moment
%           generating function.
%
%INPUTS: lambda The rate parameter of the distribution. lambda>0.
%     numDerivs The number of derivatives to take with respect to the
%               argument of the cumulant generating function, numDerivs>=0.
%             t The numPointsX1 or 1XnumPoints vector of points where the
%               moment generating function should be evaluated. If this
%               parameter is omitted or an empty matrix is passed, the
%               default value of 0 is used. This function is only valid
%               for t<lambda.
%
%OUTPUTS: cumVal A numPointsX1 or 1XnumPoints vector of the values of the
%                derivatives of the cumulant generating function given at
%                the points in t or at t=0 if t is omitted.
%
%The moment generating function of a random variable is defined to be
%E(exp(t*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter. The cumulant generating function is
%defined as the natural logarithm of the moment generating function. It can
%be shown that the moment generating function of the exponential
%distribution is
%E(exp(t*x))=lambda/(lambda-t)
%Thus, the cumulant generating function is just
%log(lambda)-log(lambda-t)
%The first derivative is thus 1/(lambda-t) and it is not dificult to deriva
%a general expression fro the nth derivative.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(t))
       t=0; 
    end

    if(numDerivs==0)
        cumVal=log(lambda)-log(lambda-t);
    else
        cumVal=factorial(numDerivs-1)/(lambda-t)^numDerivs;
    end
end

function vals=rand(N,lambda)
%%RAND Generate exponentially distributed random variables with the
%      specified rate parameter.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand 
%          returns an MXN1 matrix of random variables.
%   lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated exponential random variables.
%
%This is an implementation of the inverse transform algorithm of Chapter
%5.1 of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    U=rand(dims);

    vals=ExponentialD.invCDF(U,lambda);
end
    
function entropyVal=entropy(lambda)
%%ENTROPY Obtain the differential entropy of the exponential distribution
%         given in nats. The differential entropy of a continuous
%         distribution is entropy=-int_x p(x)*log(p(x)) dx where the
%         integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: lambda The rate parameter of the distribution. lambda>0.
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
    
    entropyVal=1-log(lambda);
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
