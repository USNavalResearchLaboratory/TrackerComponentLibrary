classdef HypergeometricD
%%HYPERGEOMETRICD Functions to handle the scalar hypergeometric
%                 distribution. Given an urn containg N1+N2 objects, N2 of
%                 type 1 and N2 of type 2, choose N balls. The number of
%                 balls of type 1 that is chosen is given by a
%                 hypergeometric distribution with parameters n, N1 and N2.
%Implemented methods are: mean, var, PMF, CDF, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release. 
    
methods(Static)
function val=mean(n,N1,N2)
%%MEAN Obtain the mean of the hypergeometric distribution.
% 
%INPUTS: n The total number of objects drawn. n<=N1+N2.
%       N1 The number of objects of type 1.
%       N2 The number of objectes of type 2. 
%
%OUTPUTS: val The mean of the hypergeometric distribution under
%             consideration.   
%
%The hypergeometric distribution is given in Chapter 2.8 of [1].
%
%EXAMPLE:
%We verify that the computed variance by comparing it to the sample
%variance.
% n=20;
% N1=63;
% N2=28;
% varVal=HypergeometricD.mean(n,N1,N2)
% numSamp=1e3;
% sampVarVal=mean(HypergeometricD.rand([1,numSamp],n,N1,N2))
%One will find both values are about 13.836
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=n*N1/(N1+N2);
end
    
function val=var(n,N1,N2)
%%VAR Obtain the variance of the hypergeometric distribution.
% 
%INPUTS: n The total number of objects drawn. n<=N1+N2.
%       N1 The number of objects of type 1.
%       N2 The number of objectes of type 2. 
%
%OUTPUTS: val The variance of the hypergeometric distribution under
%             consideration.   
%
%The hypergeometric distribution is given in Chapter 2.8 of [1].
%
%EXAMPLE:
%We verify that the computed variance by comparing it to the sample
%variance.
% n=20;
% N1=63;
% N2=28;
% varVal=HypergeometricD.var(n,N1,N2)
% numSamp=1e3;
% sampVarVal=var(HypergeometricD.rand([1,numSamp],n,N1,N2))
%One will find both values are about 3.36.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=n*N1*N2/(N1+N2)^2*(1-(n-1)/(N1+N2-1));
end
    
function val=PMF(x,n,N1,N2)
%%PMF Evaluate the hypergeometric probability mass function (PMF) at given
%     points.
%
%INPUTS: x The point(s) at which the hypergeometric PMF is evaluated.
%        n The total number of objects drawn. n<=N1+N2.
%        N1 The number of objects of type 1.
%        N2 The number of objectes of type 2.
%
%OUTPUTS: val The PMF of the hypergeometric distribution evaluated at the
%             desired point(s).
%
%The hypergeometric distribution is given in Chapter 2.8 of [1].
%
%EXAMPLE:
%Here, we validate the PMF by generating random samples and comparing the
%PMF plot with a histogram of the random samples.
% n=20;
% N1=25;
% N2=15;
% numSamples=1e5;
% figure(1)
% clf
% histogram(HypergeometricD.rand([numSamples,1],n,N1,N2),'BinWidth',1,'Normalization','pdf')
% hold on
% x=0:20;
% vals=HypergeometricD.PMF(x,n,N1,N2);
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
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(x));
    numEls=numel(x);

    denom=binomial(N1+N2,n);
    for curEl=1:numEls
        i=x(curEl);
        if(i>0&&i<=n)
            val(curEl)=binomial(N1,i).*binomial(N2,n-i)/denom;
        end
    end
end
    
function val=CDF(x,n,N1,N2)
%%CDF Evaluate the cumulative distribution function of the hypergeometric
%     distribution at desired points.
%
%INPUTS: x The point(s) at which the hypergeometric CDF is evaluated.
%        n The total number of objects drawn. n<=N1+N2.
%        N1 The number of objects of type 1.
%        N2 The number of objectes of type 2.
%
%OUTPUTS: val The CDF of the hypergeometric distribution evaluated at the
%             desired point(s).
%
%The hypergeometric distribution is given in Chapter 2.8 of [1].
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% n=20;
% N1=25;
% N2=15;
% x=12;
% numSamples=1e5;
% prob=HypergeometricD.CDF(x,n,N1,N2)
% probSamp=mean(HypergeometricD.rand([numSamples,1],n,N1,N2)<=x)
%One will find the values of both be about 0.5.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=zeros(size(x));
    numEls=numel(x);

    denom=binomial(N1+N2,n);
    for curEl=1:numEls
        i=floor(x(curEl));
        if(i>0&&i<=n)
            %Just sum the numerators in the PDF and normalize. This is
            %slow if N1+N2 is large.

            sumVal=0;
            for k=0:i
                sumVal=sumVal+binomial(N1,k).*binomial(N2,n-k);
            end
            val(curEl)=sumVal/denom;
        elseif(i>n)
           val(curEl)=1; 
        end
    end
end
    
function x=rand(N,n,N1,N2)
%%RAND Generate hypergeometric random variables.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%        n The total number of objects drawn. n<=N1+N2.
%        N1 The number of objects of type 1.
%        N2 The number of objectes of type 2.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated hypergeometric random variables.
%
%The inverse transform method fo Chapter 4.1 of [1] is used.
%
%REFERENCES:
%[1] S. M. Ross, Simulation, 4th ed. Amsterdam: Academic Press, 2006.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(isscalar(N))
        dims=[N N];
    else
        dims=N;
    end

    denom=binomial(N1+N2,n);

    x=zeros(dims);
    numEls=prod(dims);
    
    for curEl=1:numEls
        U=rand();
        cumProb=0;
        for k=0:n
            cumProb=cumProb+binomial(N1,k).*binomial(N2,n-k)/denom;

            if(U<=cumProb)
                break;
            end
        end
        x(curEl)=k;
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
