classdef BetaD
%%BETAD Functions to handle the beta distribution of the first kind.
%Implemented methods are: mean, var, PDF, CDF, rand, entropy
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
methods(Static)
function val=mean(a,b)
%%MEAN Obtain the mean of the beta distribution for given shape parameters.
%
%INPUTS: a The shape parameter that is the exponent of x; a>0.
%        b The shape parameter that is the exponent of 1-x; b>0.
%
%OUTPUTS: val The mean of the beta distribution.
%
%The mean of the beta distribution is given on the inside of the front
%cover of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=a/(a+b);
end

function val=var(a,b)
%%VAR Obtain the variance of the beta distribution for given shape
%    parameters.
%   
%INPUTS: a The shape parameter that is the exponent of x; a>0.
%        b The shape parameter that is the exponent of 1-x; b>0.
%
%OUTPUTS: val The variance of the beta distribution.
%
%The variance of the beta distribution is given on the inside of the front
%cover of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=a*b/((a+b)^2*(a+b+1));
end

function val=PDF(x,a,b)
%%PDF Evaluate the probability density function (PDF) of the beta
%     distribution at one or more points.
%
%INPUTS: x The point(s) at which the beta PDF is to be evaluated. The
%          support of the beta distribution is (0,1).
%        a The shape parameter that is the exponent of x; a>0.
%        b The shape parameter that is the exponent of 1-x; b>0. 
%
%OUTPUTS: val The value(s) of the beta PDF.
%
%The PDF of the beta distribution is gien on the inside of the front
%cover of [1].
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    %This is equivalent to val=(1/beta(a,b))*x.^(a-1).*(1-x).^(b-1);
    val=exp(-betaln(a,b)+(a-1)*log(x)+(b-1)*log(1-x));

    val(x<0|x>1)=0;
end


function val=CDF(x,a,b)
%%CDF Evaluate the cumulative distribution function (PDF) of the beta
%     distribution at one or more points.
%
%INPUTS: x The point(s) at which the beta CDF is to be evaluated. The
%          support of the beta distribution is (0,1).
%        a The shape parameter that is the exponent of x; a>0.
%        b The shape parameter that is the exponent of 1-x; b>0. 
%
%OUTPUTS: val The value(s) of the beta CDF.
%
%The PDF of the beta distribution is gien on the inside of the front
%cover of [1]. The CDF of the beta distribution is the integral from 0 to x
%of the PDF. This is just the function betainc in Matlab.
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

    selZero=(x<=0);
    selOne=(x>=1);
    val=zeros(size(x));
    val(selOne)=1;
    val(selZero)=0;

    selRest=~selZero&~selOne;   
    val(selRest)=betainc(x(selRest),a,b);
end

function val=rand(N,a,b)
%%RAND Generate beta distributed random variables with the given
%      parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%        a The shape parameter that is the exponent of x; a>0.
%        b The shape parameter that is the exponent of 1-x; b>0.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated beta random variables.
%
%The relationship between beta and gamma random variables is given in
%Example 6-12 of Chapter 6.2 of [1]. That relation is used here to generate
%the beta random variables.
%
%EXAMPLE:
%Here, we see that a histogram of the samples matches well with the PDF. 
% N=1e4;
% numPlotPoints=500;
% a=5;
% b=7;
% 
% ySamp=BetaD.rand([N,1],a,b);
% 
% x=linspace(0,1,numPlotPoints);
% y=BetaD.PDF(x,a,b);
% 
% figure(1)
% clf
% hold on
% histogram(ySamp,'Normalization','pdf')
% plot(x,y,'linewidth',2);
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%
%REFERENCES:
%[1] A. Papoulis and S. U. Pillai, Probability, Random Variables and
%    Stochastic Processes, 4th ed. Boston: McGraw Hill, 2002.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end
    
    X=GammaD.rand(dims,a,1);
    Y=GammaD.rand(dims,b,1);
    
    val=X./(X+Y);
    %It is unlikely that X and Y will both be zero, but the line below
    %corrects for the just in case.
    val(isnan(val))=0;
end

function entropyVal=entropy(a,b)
%%ENTROPY Obtain the differential entropy of the Beta distribution given in
%         nats. The differential entropy of a continuous distribution is
%         entropy=-int_x p(x)*log(p(x))  dx where the integral is over all
%         values of x. Units of nats mean that the natural logarithm is
%         used in the definition. Unlike the Shannon entropy for discrete
%         variables, the differential entropy of continuous variables can
%         be both positive and negative.
%
%INPUTS: a The shape parameter that is the exponent of x; a>0.
%        b The shape parameter that is the exponent of 1-x; b>0. 
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
    
    entropyVal=betaln(a,b)-(a-1)*psi(a)-(b-1)*psi(b)+(a+b-2)*psi(a+b);
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
