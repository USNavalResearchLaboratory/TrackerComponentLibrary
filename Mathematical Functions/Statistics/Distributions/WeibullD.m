classdef WeibullD
%Functions to handle the Weibull distribution.
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=mean(lambda,k)
%%MEAN Obtain the mean of the Weibull distirbution for given scale and
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
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

   val=lambda*gamma(1+1/k); 
end

function val=var(lambda,k)
%%VAR Obtain the variance of the Weibull distirbution for given scale and
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
%INPUTS:    x   The point(s) at which the Weibull PDF is to be evaluated.
%        lambda The scale parameter of the distribution. lambda>0.
%             k The shape parameter of the distribution. k>0.
%
%OUTPUTS:   val The value(s) of the Weibull PDF.
%
%The PDF of the Weibull distribution is given on the inside of the front
%cover of [1].
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
%     distirbution at one or more desired points.
%
%INPUTS:    x   The point(s) at which the Weibull CDF is to be evaluated.
%        lambda The scale parameter of the distribution. lambda>0.
%             k The shape parameter of the distribution. k>0.
%
%OUTPUTS:   val The value(s) of the Weibull CDF.
%
%The CDF is the integral from 0 to x of the PDF. The Weibull PDF is
%fairely simple and the integral is not hard. The Weibull PDF is given on
%the inside of the front cover of [1].
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
%INPUTS:   prob The probability or probabilities (0<=prob<=1) at which the 
%               argument of the CDF is desired.
%        lambda The scale parameter of the distribution. lambda>0.
%             k The shape parameter of the distribution. k>0.
%
%OUTPUTS:   val  The argument(s) of the CDF that would give the probability
%                or probabilities in prob.
%
%The CDF of the Weibull distribution can be easily algebraicly inverted,
%which is what is done here. The Weibull CDF is the integral from 0 to x of
%the Weibull PDF, which is given on the inside of the front cover of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=lambda*(-log(1-prob)).^(1/k);

end
    

function vals=rand(N,lambda,k)
%%RAND Generate Weibull distributed random variables with the given
%      parameters.
%
%INPUTS:    N   If N is a scalar, then rand returns an NXN matrix of random
%               variables. If N=[M,N1] is a two-element row vector, then
%               rand returns an MXN1 matrix of random variables.
%        lambda The scale parameter of the distribution. lambda>0.
%             k The shape parameter of the distribution. k>0.
%
%OUTPUTS:   vals   A matrix whose dimensions are determined by N of the
%                  generated Weibull random variables.
%
%This is an implementation of the inverse transform algorithm of Chapter
%5.1 of [1].
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

    U=rand(dims);
    vals=WeibullD.invCDF(U,lambda,k);
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
