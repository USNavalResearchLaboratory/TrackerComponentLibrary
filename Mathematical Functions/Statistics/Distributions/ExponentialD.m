classdef ExponentialD
%Functions to handle the exponential distribution.
%Implemented methods are: mean, var, PDF, CDF, invCDF, momentGenFun,
%                         cumGenFun, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(lambda)
%%VAR   Obtain the mean of the exponential distribution for a given
%       rate parameter.
%
%INPUTS:    lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS: val  The mean of the exponential distribution.
%
%The mean of the exponential distribution is given in Chapter 2.9 of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.    
    
    val=1/lambda;
end

function val=var(lambda)
%%VAR   Obtain the variance of the exponential distribution for a given
%       rate parameter.
%
%INPUTS:    lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS: val  The variance of the exponential distribution.
%
%The variance of the exponential distribution is given in Chapter 2.9 of
%[1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.    
    
    val=1/lambda^2;
end

function val=PDF(x,lambda)
%%PDF          Evaluate the PDF of the exponential probability distribution
%              function at one or more desired points.
%
%INPUTS:    x   The point(s) at which the expoenntial distribution is to be 
%               evaluated.
%        lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS:   val The value(s) of the exponential PDF with given rate
%                parameter evaluated at x.
%
%The PDF of the exponential distribution is given in Chapter 2.9 of [1].
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
%%CDF       Evaluate the CDF of the Nakagami-m probability distribution
%           function at one or more desired points.
%
%INPUTS:    x   The point(s) at which the Nakagami-m CDF is to be 
%               evaluated. Note that x>=0.
%        lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS:   prob The value(s) of the expoenntial CDF with given rate
%                parameter evaluated at x.
%
%The CDF of the exponential distribution is given in Chapter 2.9 of [1].
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
%INPUTS:    prob The probability or probabilities (0<=prob<=1) at which the 
%                argument of the CDF is desired.
%         lambda The rate parameter of the distribution. lambda>0.
%
%OUTPUTS:   val  The argument(s) of the CDF that would give the probability
%                or probabilities in prob.
%
%The CDF of the exponential distribution is very simple and is easily
%inverted using logarithms, as is done here. The CDF of the exponential
%distribution is given in Chapter 2.9 of [1].
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
%INPUTS: lambda  The rate parameter of the distribution. lambda>0.
%     numDerivs  The number of derivatives to take with respect to the
%                argument of the cumulant generating function.
%                numDerivs>=0.
%             t  The numPointsX1 or 1XnumPoints vector of points where the
%                moment generating function should be evaluated. If this
%                parameter is omitted or an empty matrix is passed, the
%                default value of 0 is used. This function is only valid
%                for t<lambda.
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
%%RAND          Generate exponentially distributed random variables with 
%               the specified rate parameter.
%
%INPUTS:    N      If N is a scalar, then rand returns an NXN 
%                  matrix of random variables. If N=[M,N1] is a two-element  
%                  row vector, then rand  returns an MXN1 matrix of 
%                  random variables.
%          lambda  The rate parameter of the distribution. lambda>0.
%
%OUTPUTS:   vals   A matrix whose dimensions are determined by N of the
%                  generated exponential random variables.
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
