classdef NegativeBinomialD
%%NEGATIVEBINOMIALD Functions to handle the negative binomial
%           distribution. The Polya distribution and the Pascal
%           distribution are special cases of this distribution. The
%           negative binomial distribution is the number of sucesses in
%           independent trials before a specified number of failures
%           occurs.
%Implemented methods are: mean, var, PMF, CDF, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
methods(Static)

function val=mean(r,p)
%%MEAN Obtain the mean of the negative binomial distribution.
%    
%INPUTS: r This is the integer number of failed trials until sampling stops
%          (defining the distribution).
%        p This is the probability of success for each trial.
%
%OUTPUTS: val The mean of the distribution.
%
%The negative binomial distribution is discussed in Chapter 4 of [1].
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% p=0.8;
% r=4;
% meanVal=NegativeBinomialD.mean(r,p)
% numSamp=1e5;
% sampMeanVal=mean(NegativeBinomialD.rand([1,numSamp],r,p))
%One will find both values are about 16.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=p*r/(1-p);
end

function val=var(r,p)
%%VAR Obtain the variance of the negative binomial distribution.
%    
%INPUTS: r This is the integer number of failed trials until sampling stops
%          (defining the distribution).
%        p This is the probability of success for each trial.
%
%OUTPUTS: val The variance of the distribution.
%
%The negative binomial distribution is discussed in Chapter 4 of [1].
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% p=0.8;
% r=4;
% varVal=NegativeBinomialD.var(r,p)
% numSamp=1e5;
% sampVarVal=var(NegativeBinomialD.rand([1,numSamp],r,p))
%One will find both values are about 80.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=p*r/(1-p)^2;
end

function val=PMF(k,r,p)
%%PMF Evaluate the probability mass function (PMF) of the negative binomial
%     distribution at given points.
%
%INPUTS: x The point(s) at which the negative binomial PMF is to be
%          evaluated. x is an integer. For nonzero PMF values, x>=r.
%        r This is the integer number of failed trials until sampling stops
%          (defining the distribution).
%        p This is the probability of success for each trial.
%
%OUTPUTS: val The value(s) of the negative binomial PMF.
%
%The negative binomial distribution is discussed in Chapter 4 of [1].

%EXAMPLE:
%Here, we validate the PMF by generating random samples and comparing the
%PMF plot with a histogram of the random samples.
% p=0.8;
% r=4;
% numSamples=1e5;
% figure(1)
% clf
% histogram(NegativeBinomialD.rand([numSamples,1],r,p),'BinWidth',1,'Normalization','pdf')
% hold on
% x=0:80;
% vals=NegativeBinomialD.PMF(x,r,p);
% stem(x,vals,'linewidth',2)
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
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(k));
    numEls=numel(val);
    
    constVal=(1-p)^r;
    for curEl=1:numEls
        if(k(curEl)>=0)
            val(curEl)=binomial(k(curEl)+r-1,k(curEl))*constVal*p^k(curEl);
        end
    end
end

function val=CDF(k,r,p)
%%CDF Evaluate the cumulative distribution function (CDF) of the negative
%     binomial distribution at given points.
%
%INPUTS: x The point(s) at which the negative binomial CDF is to be
%          evaluated. x is an integer. For nonzero CDF values, x>=r.
%        r This is the integer number of failed trials until sampling stops
%          (defining the distribution).
%        p This is the probability of success for each trial.
%
%OUTPUTS: val The value(s) of the negative binomial CDF.
%
%The negative binomial distribution is discussed in Chapter 4 of [1]. The
%CDF can be expressed in terms of a regularized incomplete beta function.
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% p=0.8;
% r=5;
% x=16;
% numSamples=1e5;
% prob=NegativeBinomialD.CDF(x,r,p)
% probSamp=mean(NegativeBinomialD.rand([numSamples,1],r,p)<=x)
%One will find the values ot both be about 0.414.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 
    
    val=1-betainc(p,k+1,r);
end

function vals=rand(N,r,p)
%%RAND Generate negative binomial random variables.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%        r This is the integer number of failed trials until sampling
%          stops (defining the distribution); r>0.
%        p This is the probability of success for each trial.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated negative binomial random variables.
%
%The negative binomial distribution arises from a Poisson distribution
%whose mean is a random variable having a central gamma distribution. This
%relationship is briefly mentioned in [1] and can be easily derived. Thus,
%this function generates gamma random variables and then generates poisson
%random variables.
%
%REFERENCES:
%[1] M. H. Quenouille, "A relation between the logarithmic, Poisson, and
%    negative binomial distributions," Biometrics, vol. 5, no. 2, pp.
%    162-164, Jun. 1949.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    lambda=GammaD.rand(dims,r,p/(1-p));

    numEls=numel(lambda);
    
    vals=zeros(size(lambda));
    for curEl=1:numEls
        vals(curEl)=PoissonD.rand(1,lambda(curEl));
    end
end
end
end