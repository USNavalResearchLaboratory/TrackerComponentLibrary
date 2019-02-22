classdef LogisticD
%%LOGISTICD Functions to handle the logistic distribution. The logistic
%    distribution is similar to the scalar normal distribution, but has
%    heavier tails.
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=mean(mu)
%%MEAN Obtain the mean of the scalar logistic probability distribution for
%      the given mean parameter.
%
%INPUTS: mu The mean (location parameter) of the distribution.
%
%OUTPUTS: val The mean of the distribution. This just returns mu, since the
%             distribution is parameterized by its mean.
%
%The logistic distribution is described in Chapter 22 of [1].
%
%REFERENCES:
%[1] N. Balakrishnan and V. B. Nevzorov, "A primer on statistical
%    distributions," 2003.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=mu;
end

function val=var(s)
%%VAR Obtain the variance of the scalar logistic probability distribution
%     for the given scale parameter.
%
%INPUTS: s The scale parameter of the distribution s>=0.
%
%OUTPUTS: val The variance of the logistic distribution under
%             consideration.
%
%The logistic distribution is described in Chapter 22 of [1]. Note that the
%definition of the scale parameter used here is the more standard
%definition so that the variance is s^2*pi^2/3, rather than the alternative
%definition where s^2 itself is the variance.
%
%REFERENCES:
%[1] N. Balakrishnan and V. B. Nevzorov, "A primer on statistical
%    distributions," 2003.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=s^2*pi^2/3;
end


function val=PDF(x,mu,s)
%%PDF Evaluate the scalar logistic probability distribution function (PDF)
%     at one or more desired points.
%
%INPUTS: x The point(s) at which the logistic PDF is to be evaluated.
%       mu The mean (location parameter) of the distribution.
%        s The scale parameter of the distribution s>=0.
%
%OUTPUTS: val The value(s) of logistic PDF evaluated at x.
%
%The logistic distribution is described in Chapter 22 of [1]. Note that the
%definition of the scale parameter used here is the more standard
%definition so that the variance is s^2*pi^2/3, rather than the alternative
%definition where s^2 itself is the variance.
%
%REFERENCES:
%[1] N. Balakrishnan and V. B. Nevzorov, "A primer on statistical
%    distributions," 2003.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=1/(4*s)*sech((x-mu)./(2*s)).^2;
end


function prob=CDF(x,mu,s)
%%CDF Evaluate the cumulative distribution function (CDF) of the scalar
%     logistic distribution at desired points.
%
%INPUTS: x The point(s) at which the logistic PDF is to be  evaluated.
%       mu The mean (location parameter) of the distribution.
%        s The scale parameter of the distribution s>=0.
%
%OUTPUTS: prob The value(s) of the CDF of the logistic distribution
%              evaluated at x.
%
%The logistic distribution is described in Chapter 22 of [1]. Note that the
%definition of the scale parameter used here is the more standard
%definition so that the variance is s^2*pi^2/3, rather than the alternative
%definition where s^2 itself is the variance.
%
%REFERENCES:
%[1] N. Balakrishnan and V. B. Nevzorov, "A primer on statistical
%    distributions," 2003.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    prob=(1/2)*(1+tanh((x-mu)./(2*s)));
end

function val=invCDF(prob,mu,s)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the scalar logistic distribution.
%    
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%          mu The mean (location parameter) of the distribution.
%           s The scale parameter of the distribution s>=0.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The logistic distribution is described in Chapter 22 of [1]. Note that the
%definition of the scale parameter used here is the more standard
%definition so that the variance is s^2*pi^2/3, rather than the alternative
%definition where s^2 itself is the variance.
%
%The CDF of the logistic distribution is simple, so the inverse CDF can be
%found by just inverting the equation.
%
%REFERENCES:
%[1] N. Balakrishnan and V. B. Nevzorov, "A primer on statistical
%    distributions," 2003.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=mu-2*s*atanh(1-2*prob);
end

function vals=rand(N,mu,s)
%%RAND Generate logistic distributed random variables with the given
%      parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element  row vector, then rand
%          returns an MXN1 matrix of random variables.
%       mu The mean (location parameter) of the distribution.
%        s The scale parameter of the distribution s>=0.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated logistic random variables.
%
%The algorithm is an implementation of the inverse transform algorithm of
%Chapter 5.1 of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end
    
    U=rand(dims);
    vals=LogisticD.invCDF(U,mu,s);
end

function entropyVal=entropy(s)
%%ENTROPY Obtain the differential entropy of the scalar logistic
%         distribution given in nats. The differential entropy of a
%         continuous distribution is entropy=-int_x p(x)*log(p(x)) dx where
%         the integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: s The scale parameter of the distribution s>=0.
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

    entropyVal=log(s)+2;
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
