classdef ParetoTypeID
%%PARETOTYPID Functions to handle type I of the Pareto distribution.
%Implemented methods are: PDF, CDF, invCDF, mean, var, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=PDF(x,xMin,a)
%%PDF Evaluate the probability density function (PDF) of the type I Pareto
%     distribution.
%
%INPUTS: x The matrix of points at which the PDF should be evaluated.
%     xMin The scale parameter of the distribution. This marks the lowest
%          point having nonzero likelihood.
%        a The shape parameter of the distribution. a>0.
%
%OUTPUTS: val The value of the PDF evaluated at the point(s) in x.
%
%The PDF of the distribution is given in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Pareto Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/ParetoDistribution.html
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=a*xMin^a./x.^(a+1);
    val(x<xMin)=0;
end

function val=CDF(x,xMin,a)
%%CDF Evaluate the cumulative distribution function (CDF) of the type I
%     Pareto distribution.
%
%INPUTS: x The matrix of points at which the CDF should be evaluated.
%     xMin The scale parameter of the distribution. This marks the lowest
%          point having nonzero likelihood.
%        a The shape parameter of the distribution. a>0.
%
%OUTPUTS: val The value of the CDF evaluated at the point(s) in x.
%
%The CDF of the distribution is given in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Pareto Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/ParetoDistribution.html
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=1-(xMin./x).^a;
    val(x<xMin)=0;
end

function val=invCDF(prob,xMin,a)
%%INVCDF          Evaluate the inverse of the cumulative distribution
%                 function of type I of the Pareto distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%        xMin The scale parameter of the distribution. This marks the
%             lowest point having nonzero likelihood.
%           a The shape parameter of the distribution. a>0.
%
%OUTPUTS:   val  The argument(s) of the CDF that would give the probability
%                or probabilities in prob.
%
%The CDF in [1] is simple and the algebraic inverse can be easily found.
%Note that since values below xMin all have 0 probability, the inverse of 0
%is assigned just to xMin.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Pareto Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/ParetoDistribution.html
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

val=(1-prob).^(-1/a)*xMin;
end

function val=mean(xMin,a)
%%MEAN Evaluate the mean of the type I Pareto distribution.
%
%INPUTS: xMin The scale parameter of the distribution. This marks the
%             lowest point having nonzero likelihood.
%           a The shape parameter of the distribution. a>0.
%
%OUTPUTS: val The mean of the distribution with the given parameters.
%
%The mean of the distribution is given in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Pareto Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/ParetoDistribution.html
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(a<=1)
        val=Inf;
    else
        val=a*xMin/(a-1);
    end
end
    
function val=var(xMin,a)
%%MEAN Evaluate the variance of the type I Pareto distribution.
%
%INPUTS: xMin The scale parameter of the distribution. This marks the
%             lowest point having nonzero likelihood.
%           a The shape parameter of the distribution. a>0.
%
%OUTPUTS: val The mean of the distribution with the given parameters.
%
%The variance of the distribution is given in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Pareto Distribution." From MathWorld--A Wolfram
%    Web Resource. http://mathworld.wolfram.com/ParetoDistribution.html
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(a<=2)
       val=Inf; 
    else
        val=a*xMin^2/((a-1)^2*(a-2));
    end
end

function vals=rand(N,xMin,a)
%%RAND Generate Pareto type-I distributed random variables with the given
%      parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%     xMin The scale parameter of the distribution. This marks the
%          lowest point having nonzero likelihood.
%        a The shape parameter of the distribution. a>0.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Pareto type I random variables.
%
%The algorithm is an implementation of the inverse transform algorithm of
%Chapter 5.1 of [1]. When the noncentral distribution is used, the random
%variables are generated by summing the squares of normally distributed
%random variables.
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
    vals=ParetoTypeID.invCDF(U,xMin,a);
end
    
function entropyVal=entropy(xMin,a)
%%ENTROPY Obtain the differential entropy of the Pareto type I distribution
%         given in nats. The differential entropy of a continuous
%         distribution is entropy=-int_x p(x)*log(p(x)) dx where the
%         integral is over all values of x. Units of nats mean that the
%         natural logarithm is used in the definition. Unlike the Shannon
%         entropy for discrete variables, the differential entropy of
%         continuous variables can be both positive and negative.
%
%INPUTS: xMin The scale parameter of the distribution. This marks the
%             lowest point having nonzero likelihood.
%           a The shape parameter of the distribution. a>0.
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

    entropyVal=log(xMin/a)+1+1/a;
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
