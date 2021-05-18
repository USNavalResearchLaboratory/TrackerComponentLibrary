classdef PoissonD
%%POISSOND Functions to handle the Poisson distribution.
%Implemented methods are: mean, var, PMF, CDF, momentGenFun, cumGenFun,
%                         rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(lambda)
%%MEAN Obtain the mean of the Poisson distribution for the specified mean
%      parameter.
%
%INPUTS: lambda The mean (and variance) of the Poisson distribution.
%
%OUTPUTS: val The mean of the Poisson distribution under consideration.
%
%The Poisson distribution is discussed in Chapter 2.8 of [1].
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% lambda=10;
% meanVal=PoissonD.mean(lambda)
% numSamp=1e5;
% sampMeanVal=mean(PoissonD.rand([1,numSamp],lambda))
%One will find both values are about 10.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=lambda;
end

function val=var(lambda)
%%VAR Obtain the variance of the Poisson distribution for the specified
%     mean parameter
%
%INPUTS: lambda The mean (and variance) of the Poisson distribution,
%
%OUTPUTS: val The variance of the Poisson distribution.
%
%The Poisson distribution is discussed in Chapter 2.8 of [1].
%
%EXAMPLE:
%We verify that the computed variance by comparing it to the sample
%variance.
% lambda=10;
% varVal=PoissonD.var(lambda)
% numSamp=1e5;
% sampVarVal=var(PoissonD.rand([1,numSamp],lambda))
%One will find both values are about 10.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=lambda;
end

function val=PMF(x,lambda)
%%PMF Evaluate the Poisson probability mass function (PMF) at given points
%     with a specified mean parameter.
%
%INPUTS: x The point(s) at which the Poisson PMF is to be evaluated. x is
%          an integer. For nonzero PMF values, x>=0.
%   lambda The mean (and variance) of the Poisson distribution under
%          consideration.
%
%OUTPUTS: val The value(s) of the Poisson PMF with mean lambda.
%
%The Poisson distribution is discussed in Chapter 2.8 of [1].
%
%EXAMPLE:
%Here, we validate the PMF by generating random samples and comparing the
%PMF plot with a histogram of the random samples.
% lambda=6;
% numSamples=1e5;
% figure(1)
% clf
% histogram(PoissonD.rand([numSamples,1],lambda),'BinWidth',1,'Normalization','pdf')
% hold on
% x=0:16;
% vals=PoissonD.PMF(x,lambda);
% stem(x,vals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    val=zeros(size(x));
    sel=(x>=0);

    val(sel)=(lambda.^x(sel)./factorial(x(sel)))*exp(-lambda);
end

function [FVal,pVal]=CDF(x,lambda)
%%CDF Evaluate the cumulative distribution function of the Poisson
%     distribution at desired points.
%
%INPUTS: x The integer point(s) at which the Poisson CDF is evaluated.
%   lambda The mean (and variance) of the Poisson distribution.
%
%OUTPUTS: FVal The CDF of the Poisson distribution evaluated at the desired
%              point(s).
%         pVal The value of the PMF of the Poisson distribution at the
%              desired point(s).
%
%The recursion is from Chapter 4.2 of [1].
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% lambda=6;
% x=3;
% numSamples=1e5;
% prob=PoissonD.CDF(x,lambda)
% probSamp=mean(PoissonD.rand([numSamples,1],lambda)<=x)
%One will find the values ot both be about 0.151.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    FVal=zeros(size(x));
    pVal=zeros(size(x));
    
    x=floor(x);%In case decimals were passed.
    numEls=numel(x);

    for curEl=1:numEls
        p=exp(-lambda);
        F=p;
        for I=1:fix(x(curEl))
            p=p*lambda/I;
            F=F+p;
        end
        FVal(curEl)=F;
        pVal(curEl)=p;
    end
    FVal(x<0)=0;
    pVal(x<0)=0;
end

function momentVal=momentGenFun(lambda,numDerivs,t)
%%MOMENTGENFUN Evaluate the moment generating function (or one of its
%              derivatives) of the Poisson distribution. Taking the kth
%              derivative of the moment generating function and evaluating
%              it at t=0 provides the kth noncentral moment of the
%              distribution.
%
%INPUTS: lambda The mean (and variance) of the Poisson distribution.
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
%generating function of the Poisson distribution is
%f^0(t)=exp(lambda*(exp(t)-1))
%By noting the pattern in differentiation (using the chain rule), it can be
%seen that the derivatives of the moment generating function can be
%recurively expressed. Specifically, if f^n(t) is the nth derivative, the
%value for n>0 is 
%f^n(t)=lambda*exp(t)*(sum_{k=0}^(n-1) binomial(n-1,k)*f^k(t))
%This function evaluates the moments using the above recursion. The
%binomial terms are computed recurively using the recursion inherent to
%Pascal's triangle.
%
%Note that because the function recursively accumulates values for the
%derivatives, a lot of memory will be used if numDerivs is very large.
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<3||isempty(t))
    t=0; 
end

numT=length(t);
momentVal=zeros(numT,1);%Allocate Space
for curT=1:numT
    tVal=t(curT);
    
    fVals=zeros(numDerivs+1,1);
    fVals(1)=exp(lambda*(exp(tVal)-1));%The moment generating function value.

    pascalTriangleRowOld=zeros(numDerivs+1,1);
    pascalTriangleRowOld(1)=1;
    expCoeff=exp(tVal)*lambda;
    for curDeriv=1:numDerivs
        %Get the next row of Pascal's triangle
        pascalTriangleRow=zeros(numDerivs+1,1);
        pascalTriangleRow(1)=1;
        for i=2:(curDeriv-1)
            pascalTriangleRow(i)=pascalTriangleRowOld(i-1)+pascalTriangleRowOld(i);
        end
        pascalTriangleRow(curDeriv)=1;

        fVals(curDeriv+1)=expCoeff*sum(pascalTriangleRow.*fVals);
        pascalTriangleRowOld=pascalTriangleRow;
    end

    momentVal(curT)=fVals(end);
end

end

function cumVal=cumGenFun(lambda,numDerivs,t)
%%CUMGENFUN Evaluate the cumulant generating function (or one of its
%           derivatives) of the Poisson distribution. The cumulant
%           generating function is the natural logarithm of the moment
%           generating function.
%
%INPUTS: lambda The mean (and variance) of the Poisson distribution.
%     numDerivs The number of derivatives to take with respect to the
%               argument of the cumulant generating function, numDerivs>=0.
%             t The numPointsX1 or 1XnumPoints vector of points where the
%               moment generating function should be evaluated. If this
%               parameter is omitted or an empty matrix is passed, the
%               default value of 0 is used.
%
%OUTPUTS: cumVal A numPointsX1 or 1XnumPoints vector of the values of the
%                derivatives of the cumulant generating function given at
%                the points in t or at t=0 if t is omitted.
%
%The moment generating function of a random variable is defined to be
%E(exp(t*x)) where E is the expected value operator, x is the random
%variable and t is a real parameter. The cumulant generating function is
%defined as the natural logarithm of the moment generating function. It can
%be shown that the moment generating function of the Poisson distribution
%is
%E(exp(t*x))=exp(lambda*(exp(t)-1))
%Thus, the cumulant generating function is just
%lambda*(exp(t)-1)
%All derivatives of the cumulant generating function are consequently
%exp(t)*lambda.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(nargin<3||isempty(t))
    t=0; 
end

if(numDerivs==0)
    cumVal=lambda*(exp(t)-1);
else
    cumVal=exp(t)*lambda;
end
    
end

function vals=rand(N,lambda)
%%RAND Generate Poisson random variables with a given mean.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%   lambda The mean (and variance) of the Poisson distribution.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Poisson random variables.
%
%The Poisson random variables are generated according to the simple
%method given in Chapter 4.2 of [1].
%
%The algorithm can handle small and large values of lambda. However, there
%is always a possibility that the value desired might be on the tail
%end of the distribution and be very slow. Thus, the number of iterations
%(the offset from the mean) is set to the maximum of 10 and 10* the
%variance.
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

    maxIter=max(10,10*lambda);

    vals=zeros(dims);
    for curRow=1:dims(1)
        for curCol=1:dims(2)
            I=fix(lambda);
            [F,p]=PoissonD.CDF(I,lambda);
            U=rand(1);

            if(U<=F)
            %This means that the random number being generated is less than  
            %or equal to I and we must search downward.
                while(U<F&&I>-1)
                    F=F-p;
                    p=p*(I/lambda);
                    I=I-1;
                end
                I=I+1;
            else
            %This means that the random number being generated is greater 
            %than I and we must search upward starting from I+1.
                curIter=0;
                while(U>F&&curIter<maxIter)
                    p=p*lambda/(I+1);
                    F=F+p;
                    I=I+1;
                    curIter=curIter+1;
                end
            end

            vals(curRow ,curCol)=I;
        end
    end
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
