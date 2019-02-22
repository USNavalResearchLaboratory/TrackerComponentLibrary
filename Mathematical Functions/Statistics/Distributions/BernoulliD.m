classdef BernoulliD
%%BERNOULLID Functions to handle the Bernoulli distribution. This is just a
%            binary distribution having a certain probability p of being 1.
%Implemented methods are: mean, var, PMF, CDF, rand, momentGenFun, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
function val=mean(p)
%%MEAN Obtain the mean of the Bernoulli distribution.
%
%INPUTS: p The probability of success of a Bernoulli trial.
%
%OUTPUTS: val The mean of the Bernoulli distribution.
%
%As noted on the page opposite the inside front cover of [1], the mean is
%just p.
%
%EXAMPLE:
%We show that the mean matches the sample mean.
% p=0.75;
% meanVal=BernoulliD.mean(p)
% numRuns=1e6;
% meanSamp=mean(BernoulliD.rand([numRuns,1],p))
%One will see that both values are about 0.75.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 
   
    val=p;
end

function val=var(p)
%%VAR Obtain the variance of the Bernoulli distribution.
%
%INPUTS: p The probability of success of a Bernoulli trial.
%
%OUTPUTS: val The variance of the Bernoulli distribution.
%
%As noted on the page opposite the inside front cover of [1], the mean is
%just p*(1-p).
%
%EXAMPLE:
%We show that the mean matches the sample mean.
% p=0.75;
% varVal=BernoulliD.var(p)
% numRuns=1e6;
% varSamp=var(BernoulliD.rand([numRuns,1],p))
%One will see that both values are about 0.1875.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 

   val=p*(1-p); 
end

function val=PMF(x,p)
%%PMF Evaluate the Bernoulli probability mass function (PMF) at given
%     points.
%
%INPUTS: x The point(s) at which the Bernoulli PMF is to be evaluated. For
%          nonzero PMF values, x=0 or x=1.
%        p The probability of success of a Bernoulli trial.
%
%OUTPUTS: val The value(s) of the Bernoulli PMF.
%
%The Bernoulli distribution is given on the page opposite the inside front
%cover of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(x));
    val(x==0)=1-p;
    val(x==1)=p;
end
    
function val=CDF(x,p)
%%CDF Evaluate the cumulative distribution function of the Bernoulli
%     distribution at desired points.
%
%INPUTS: x The point(s) at which the Bernoulli PMF is to be evaluated. For
%          nonzero PMF values, x=0 or x=1.
%        p The probability of success of a Bernoulli trial.
%
%OUTPUTS: val The CDF of the Bernoulli distribution evaluated at the
%             desired point(s).
%
%The Bernoulli distribution is given on the page opposite the inside front
%cover of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(x));
    sel=(0<=x|x<1);
    val(sel)=1-p;
    sel=(x>=1);
    val(sel)=p;
end

function x=rand(N,p)
%%RAND Generate Bernoulli random variables with a given mean.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%        p The probability of success of a Bernoulli trial.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Bernoulli random variables.
%
%The Bernoulli random variables are generated simply using the comparison
%mentioned in Chapter 4 of [1]. [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end
    
    x=rand(dims)<=p;
end

function momentVal=momentGenFun(p,numDerivs,t)
%%MOMENTGENFUN Evaluate the moment generating function (or one of its
%              derivatives) of the Bernoulli distribution. Taking the kth
%              derivative of the moment generating function and evaluating
%              it at t=0 provides the kth noncentral moment of the
%              distribution.
%
%INPUTS: p The probability of success of a Bernoulli trial.
% numDerivs The number of derivatives to take with respect to the
%          argument of the moment generating function. numDerivs>=0.
%        t The matrix of points where the moment generating function
%          should be evaluated. If this parameter is omitted or an
%          empty matrix is passed, the default value of 0 is used.
%
%OUTPUTS: momentVal A numPointsX1 vector of the values of the derivatives
%                   of the moment generating function given at the points
%                   in t or at t=0 if t is omitted.
%
%The moment generating function of a random variable is defined to be
%E(exp(t*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter. The moment generating function for the
%Bernoulli distribution is just (1-p)+p*exp(t). Derivatives are
%straightforward.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(t))
        t=0; 
    end

    momentVal=exp(t)*p;
    sel=numDerivs==0;
    momentVal(sel)=momentVal(sel)+1-p;
end

function entropyVal=entropy(p)
%%ENTROPY Obtain the Shannon entropy of the Bernoulli distribution given in
%         nats. The Shannon entropy of a discrete distribution is
%         entropy=-sum_x Pr(x)*log(Pr(x)) where the sum is over all
%         discrete values of x. Units of nats mean that the natural
%         logarithm is used in the definition.
%
%INPUTS: p The probability of success of a Bernoulli trial.
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
    
    if(p==0||p==1)
        entropyVal=0;
        return;
    end

    q=1-p;

    entropyVal=-q*log(q)-p*log(p);
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
