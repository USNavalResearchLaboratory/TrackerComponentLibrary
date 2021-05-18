classdef GeometricD
%%GEOMETRICD Functions to handle two variants of the geometric
%            distribution. A geometric random variable with probability p
%            arises when determinign the number of independent Bernoulli
%            trials that must be performed until a success is obtained. The
%            probability of success for each trial is p.
%Implemented methods are: mean, var, PMF, CDF, invCDF, momentGenFun,
%                         cumGenFun, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(p,type)
%%MEAN Obtain the mean of the geometric distribution.
%
%INPUTS: p The probability of success of each Bernoulli trial.
%     type A parameter indicating the type of Geometric distribution.
%         Possible values are:
%         0 (The default if omitted or an empty matrix is passed) The
%           distribution is of the number of Bernoulli trials needed before
%           obtaining a success (a 1).
%         1 The distribution is of the number of failures obtained from
%           Bernoulli trials before obtaining a success.
%
%OUTPUTS: val The mean of the geometric distribution under consideration.
%
%The mean of the geometric distribution, is given in table 5-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(type))
       type=0; 
    end

    if(type==0)
        val=1/p;
    elseif(type==1)
        val=(1-p)/p;
    else
        error('Unknown type specified');
    end
end

function val=var(p,type)
%%VAR Obtain the variance of the geometric distribution
%
%INPUTS: p The probability of success of each Bernoulli trial.
%     type A parameter indicating the type of Geometric distribution.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) The
%            distribution is of the number of Bernoulli trials needed
%            before obtaining a success (a 1).
%          1 The distribution is of the number of failures obtained from
%            Bernoulli trials before obtaining a success.
%
%OUTPUTS: val The variance of the geometric distribution.
%
%The variance of the geometric distribution, is given in table 5-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(type))
       type=0; 
    end

    if(type==0||type==1)
        val=(1-p)/p^2;
    else
        error('Unknown type specified');
    end
end

function val=PMF(x,p,type)
%%PMF Evaluate the geometric probability mass function (PMF) at given
%     points.
%
%INPUTS: x The point(s) at which the Poisson PMF is to be evaluated. x is
%          an integer. For a type 0 distribution x>=1, for a type 1
%          distribution, x>=0.
%        p The probability of success of each Bernoulli trial.
%     type A parameter indicating the type of Geometric distribution.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) The
%            distribution is of the number of Bernoulli trials needed
%            before obtaining a success (a 1).
%          1 The distribution is of the number of failures obtained from
%            Bernoulli trials before obtaining a success.
%
%OUTPUTS: val The value(s) of the geometric PMF
%
%The PMF of the geometric distribution, is given in table 5-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(nargin<3||isempty(type))
       type=0; 
    end

    if(type==0)
        val=(1-p).^(x-1)*p;
    elseif(type==1)
        val=(1-p).^(x)*p;
    else
        error('Unknown type specified');
    end
end

function val=CDF(x,p,type)
%%PMF Evaluate the cumulative distribution function (CDF)of the geometric
%     distribution at desired points.
%
%INPUTS: x The point(s) at which the Poisson CDF is to be evaluated. x is
%          an integer. For a type 0 distribution x>=1, for a type 1
%          distribution, x>=0.
%        p The probability of success of each Bernoulli trial.
%     type A parameter indicating the type of Geometric distribution.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) The
%            distribution is of the number of Bernoulli trials needed
%            before obtaining a success (a 1).
%          1 The distribution is of the number of failures obtained from
%            Bernoulli trials before obtaining a success.
%
%OUTPUTS: val The value(s) of the geometric CDF
%
%The CDF of the geometric distribution can be obtained by summing the PMF
%of the distribution. The PMF, is given in Table 5-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(type))
       type=0; 
    end

    if(type==0)
        val=1-(1-p).^x;
    elseif(type==1)
        val=1-(1-p).^(x+1);
    else
        error('Unknown type specified');
    end
end


function x=invCDF(prob,p,type)
%%INVCDF Evaluate the inverse of the cumulative distribution function of
%        the geometric distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%           p The probability of success of each Bernoulli trial.
%        type A parameter indicating the type of Geometric distribution.
%             Possible values are:
%             0 (The default if omitted or an empty matrix is passed) The
%               distribution is of the number of Bernoulli trials needed
%               before obtaining a success (a 1).
%             1 The distribution is of the number of failures obtained from
%               Bernoulli trials before obtaining a success.
%
%OUTPUTS: val The argument(s) of the CDF that would give the probability or
%             probabilities in prob.
%
%The CDF of the geometric distribution is very simple and is easily
%inverted using logarithms, as is done here. However, since it is a
%discrete distribution, the continuous result is rounded up to the next
%highest integer if type=0 (and forced to be a minimum of 1) and truncated
%to the next lowest integer if type=1 (and forced to be a minimum of 0).
%The PMF of the geometric distribution, from which the CDF
%can be derived, is given in Table 5-2 of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(type))
       type=0; 
    end

    if(type==0)
        x=log(1-prob)/log(1-p);
        %If any of the results are NaN, then it should be because p=1
        %(meaning success on first trial) and prob=1. In this degenerate
        %case, x=1 --a single trial is needed until success.
        sel=isnan(x);
        x(sel)=1;
        
        %Round to next highest integer with the minimum value being 1.
        x=max(ceil(x),1);
    elseif(type==1)
        x=log(1-prob)/log(1-p)-1;
        %If any of the results are NaN, then it should be because p=1
        %(meaning success on first trial) and prob=1. In this degenerate
        %case, x=0 --no failures before success
        sel=isnan(x);
        x(sel)=0;
        
        %Round to the next highest integer with the minimum value being 0.
        x=max(fix(x),0);
    else
        error('Unknown type specified');
    end
end

function vals=rand(N,p,type)
%%RAND Generate geometrically distributed random variables.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of  random variables.
%        p The probability of success of each Bernoulli trial.
%     type A parameter indicating the type of Geometric distribution.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) The
%            distribution is of the number of Bernoulli trials needed
%            before obtaining a success (a 1).
%          1 The distribution is of the number of failures obtained from
%            Bernoulli trials before obtaining a success.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated geometric random variables.
%
%This is an implementation of the inverse transform algorithm of Chapter
%5.1 of [1]. That algorithm is for a continous distribution, but it can be
%applied here, because there is no upper limit to the value of a geometric
%random variables. It is in the same sense as the inverse transform method
%of Chapter 4.1 of [1] for discrete distributions.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(type))
       type=0; 
    end

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    U=rand(dims);

    vals=GeometricD.invCDF(U,p,type);
end

function momentVal=momentGenFun(p,type,numDerivs,t)
%%MOMENTGENFUN Evaluate the moment generating function (or one of its
%              derivatives) of the Geometric distribution. Taking the kth
%              derivative of the moment generating function and evaluating
%              it at t=0 provides the kth noncentral moment of the
%              distribution.
%
%INPUTS: p The probability of success of each Bernoulli trial.
%     type A parameter indicating the type of Geometric distribution.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) The
%            distribution is of the number of Bernoulli trials needed
%            before obtaining a success (a 1).
%          1 The distribution is of the number of failures obtained from
%            Bernoulli trials before obtaining a success.
%  numDerivs The number of derivatives to take with respect to the
%            argument of the moment generating function. numDerivs>=0.
%          t The vector or matrix of points where the moment generating
%            function should be evaluated. If this parameter is omitted or
%            an empty matrix is passed, the default value of 0 is used.
%
%OUTPUTS: momentVal A matrix of the values of the specified derivative of
%                   the moment generating function given at the points in t
%                   or at t=0 if t is omitted.
%
%The moment generating function of a random variable is defined to be
%E(exp(t*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter. It can be shown that the moment
%generating function of the geometric distribution that counts the number
%of Bernoulli trials needed to get one success is
%E(exp(t*x))=p*exp(t)/(1-(1-p)*exp(t))
%After considerable work taking derivatives, one can identify the pattern
%that the numDerivs derivative of the moment generating function is given
%by the sum
%\sum_{i=0}^n p*(1-p)^i*factorial(i)*StirlingNumber2(n+1,i+1)*(exp(t)/(1-(1-p)*exp(t)))^(i+1)
%where n is the number of derivatives to take n>=0.
%
%On the other hand, the moment generating function of the geometric
%distribution defined as the number of failed Bernoulli trials before a
%success can be shown to be
%E(exp(t*x))=p/(1-(1-p)*exp(t))
%After considerable work taking derivatives, one find the numDerivs
%derivative of the moment generating function is
%\sum_{i=0}^(numDerivs-1)p*A(n,i+1)*(1-p)^(i+1)*exp(t*(i+1))/(1-(1-p)*exp(t))^(i+2)
%where n is the number of derivatives to take, n>=1 and the function A is
%A(n,k)=\sum_{m=1}^{max(k,n)}(-1)^(m+mod(k,2))*binomial(k,m)*m^n
%Note that for k>=n, the function A(n,k) can be written as
%factorial(n)*StirlingNumber2(k,n)
%The A function is implemented below in this file.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(type))
       type=0; 
    end

    if(nargin<4||isempty(t))
       t=0; 
    end

    n=numDerivs;

    if(type==0)
        %The case where the distribution is the number of trials needed to
        %get a success. 

        sumVal=0;
        ratVal=p*exp(t)./(1-(1-p)*exp(t));
        ratMultVal=(1-p)*exp(t)./(1-(1-p)*exp(t));
        for i=0:n
            sumVal=sumVal+StirlingNumber2(n+1,i+1)*ratVal;
            %The i+1 is for the factorial term.
            ratVal=ratVal.*(i+1).*ratMultVal;
        end
        momentVal=sumVal;
        return;
    elseif(type==1)
        %The case where the distribution is the number of failed trials
        %before a success.
        if(n==0)
            momentVal=p./(1-(1-p)*exp(t));
            return;
        else
            sumVal=0;
            ratVal=p*(1-p)*exp(t)./(1-(1-p)*exp(t)).^2;
            ratMultVal=(1-p)*exp(t)./(1-(1-p)*exp(t));
            for i=0:(n-1)
                sumVal=sumVal+A(n,i+1)*ratVal;
                ratVal=ratVal.*ratMultVal;
            end
            momentVal=sumVal;
            return;
        end
    else
        error('Unknown type specified');
    end
    
end

function cumVal=cumGenFun(p,type,numDerivs,t)
%%CUMGENFUN Evaluate the cumulant generating function (or one of its
%           derivatives) of the geometric distribution. The cumulant
%           generating function is the natural logarithm of the moment
%           generating function.
%
%INPUTS: p The probability of success of each Bernoulli trial.
%     type A parameter indicating the type of Geometric distribution.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) The
%            distribution is of the number of Bernoulli trials needed
%            before obtaining a success (a 1).
%          1 The distribution is of the number of failures obtained from
%            Bernoulli trials before obtaining a success.
%  numDerivs The number of derivatives to take with respect to the
%            argument of the cumulant generating function. numDerivs>=0.
%          t The vector or matrix of points where the moment generating
%            function should be evaluated. If this parameter is omitted or
%            an empty matrix is passed, the default value of 0 is used.
%
%OUTPUTS: cumVal A  matrix of the values of the specified derivative of
%                the cumulant generating function given at the points in t
%                or at t=0 if t is omitted.
%
%The moment generating function of a random variable is defined to be
%E(exp(t*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter. The cumulant generating function is
%defined as the natural logarithm of the moment generating function. It can
%be shown that the moment generating function geometric distribution that
%counts the number of Bernoulli trials needed to get one success is
%E(exp(t*x))=p*exp(t)/(1-(1-p)*exp(t))
%Thus, the cumulant generating function is
%log(E(exp(t*x)))=log(p)+t-log(1-(1-p)*exp(t));
%After considerable work taking derivatives, one can identify the pattern
%that the numDerivs derivative of the moment generating function is given
%by the sum
%\sum_{i=0}^{n-1} (1-p)^(i+1)*factorial(i)*StirlingNumber2(n,i+1)*(exp(t)/(1-(1-p)*exp(t)))^(i+1)
%where n is the number of derivatives to take and n>=2. The first
%derivative is handled as a special case (due to the t term differentiating
%to 1).
%
%On the other hand, the moment generating function of the geometric
%distribution defined as the number of failed Bernoulli trials before a
%success can be shown to be
%E(exp(t*x))=p*exp(t)/(1-(1-p)*exp(t))
%This corresponds to a cumulant generating function of 
%log(E(exp(t*x)))=log(p)-log(1-(1-p)*exp(t))
%After considerable work taking derivatives, one find the numDerivs
%derivative of the cumulant generating function is
%\sum_{i=0}^{n}(1-p)^(i+1)*factorial(i)*StirlingNumber2(n,i+1)*(exp(t)/(1-(1-p)*exp(t)))^(i+1)
%for n>=1.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(type))
       type=0; 
    end

    if(nargin<4||isempty(t))
       t=0; 
    end

    n=numDerivs;
    
    if(type==0)
        %The case where the distribution is the number of trials needed to
        %get a success. 
        
        if(numDerivs==0)
            cumVal=log(p)+t-log(1-(1-p)*exp(t));
            return;
        elseif(numDerivs==1)
            cumVal=1+exp(t).*(1-p)./(1-(1-p)*exp(t));
            return;
        else
            sumVal=0;
            ratVal=(1-p)*exp(t)./(1-(1-p)*exp(t));
            ratMultVal=(1-p)*exp(t)./(1-(1-p)*exp(t));
            for i=0:(n-1)
                sumVal=sumVal+StirlingNumber2(n,i+1)*ratVal;
                %The i+1 is for the factorial term.
                ratVal=ratVal.*(i+1).*ratMultVal;
            end
            cumVal=sumVal;
            return;
        end
    elseif(type==1)
        %The case where the distribution is the number of failed trials
        %before a success.
        if(numDerivs==0)
            cumVal=log(p)-log(1-(1-p)*exp(t));
        else
            sumVal=0;
            ratVal=(1-p)*exp(t)./(1-(1-p)*exp(t));
            ratMultVal=(1-p)*exp(t)./(1-(1-p)*exp(t));
            for i=0:(n-1)
                sumVal=sumVal+StirlingNumber2(n,i+1)*ratVal;
                ratVal=ratVal.*(i+1).*ratMultVal;
            end

            cumVal=sumVal;
            return;
        end
    else
        error('Unknown type specified');
    end
end

function entropyVal=entropy(p,type)
%%ENTROPY Obtain the Shannon entropy of the Geometric distribution given in
%         nats. The Shannon entropy of a discrete distribution is
%         entropy=-sum_x Pr(x)*log(Pr(x)) where the sum is over all
%         discrete values of x. Units of nats mean that the natural
%         logarithm is used in the definition.
%
%INPUTS: p The probability of success of each Bernoulli trial.
%     type A parameter indicating the type of Geometric distribution.
%          Possible values are:
%          0 (The default if omitted or an empty matrix is passed) The
%            distribution is of the number of Bernoulli trials needed
%            before obtaining a success (a 1).
%          1 The distribution is of the number of failures obtained from
%            Bernoulli trials before obtaining a success. 
%
%OUTPUTS: entropyVal The value of the Shannon entropy in nats.
%
%Shannon entropy is defined in Chapter 2 of [1].
%
%REFERENCES:
%[1] T. M. Cover and J. A. Thomas, Elements of Information Theory, 2nd ed.
%    Hoboken, NJ: Wiley-Interscience, 2006.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<2||isempty(type))
       type=0; 
    end

    switch(type)
        case 0%Trials before success
            entropyVal=-((1-p)/p)*log(1-p)-log(p);
        case 1%Failures before success
            %The same as case 0.
            entropyVal=-((1-p)/p)*log(1-p)-log(p); 
        otherwise
            error('Unknown type specified.')
    end
end
end
end

function val=A(n,k)
%A This implements the function
%  A(n,k)=\sum_{m=1}^{max(k,n)}(-1)^(m+mod(k,2))*binomial(k,m)*m^n
%  which plays a role in the computation of coefficients for the moment
%  generating function. For k<n, this is not implemented in any particular
%  manner that might minimize the rish of overflow/ underflow of
%  intermediate results. For k>=n, it is equivalent to
%  factorial(n)*StirlingNumber2(k,n).

if(k>=n)
    val=factorial(n)*StirlingNumber2(k,n);
else
    val=0;
    signVal=(-1)^(1+mod(k,2));
    for m=1:max(k,n)
        val=val+signVal*binomial(k,m)*m^n;
        signVal=-signVal;
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
