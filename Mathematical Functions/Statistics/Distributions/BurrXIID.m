classdef BurrXIID
%%BURRXIID Function to handle scalar Burr type XII distributions. These are
%       continuous distributions for values >=0 that can support wide
%       ranges of skew and kurtosis values. It has found use in various
%       econometric applications. The distribution is discussed in [1].
%Implemented methods are: mean, var, PDF, CDF, invCDF, rand, kthMoment
%
%REFERENCES:
%[1] R. N. Rodriguez, "A guide to the Burr type XII distributions,"
%    Biometrica, vol. 64, no. 1, pp. 129-134, Apr. 1977.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)

function val=mean(c,k,lambda)
%%MEAN Obtain the mean of the Burke type XII distribution. If c*k>=1, then
%      the mean does not exist and an error is thrown.
%
%INPUTS: c A scalar parameter > 0 that is the exponent of x in the PDF.
%        k A scalar parameter >0 that is the exponent of the denominator in
%          the PDF-1.
%   lambda A scalar parameter >0 that x is divided by in the PDF. This is
%          an extension from the original formulation of the distribution
%          in that x can be scaled. If this parameter is omitted or an
%          empty matrix is passed, then the default value of 1 is used.
%
%OUTPUTS: val The mean of the Burke type XII distribution.
%
%An expression for arbitrary noncentral moments of the distribution is
%given in [1]. This is just a special case with a change of variables to
%deal with lambda.
%
%EXAMPLE:
%We verify the computed mean by comparing it to the sample mean.
% c=1.75;
% k=3;
% lambda=7;
% meanVal=BurrXIID.mean(c,k,lambda)
% numSamp=1e6;
% sampMeanVal=mean(BurrXIID.rand([numSamp,1],c,k,lambda))
%One will see that both values are about 3.9458.
%
%REFERENCES:
%[1] R. N. Rodriguez, "A guide to the Burr type XII distributions,"
%    Biometrica, vol. 64, no. 1, pp. 129-134, Apr. 1977.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 
    
    if(nargin<3||isempty(lambda))
       lambda=1; 
    end

    val=BurrXIID.kthMoment(1,c,k,lambda);
end

function val=var(c,k,lambda)
%%VAR Obtain the variance of the Burke type XII distribution. If c*k>=1,
%      then the variance does not exist and an error is thrown.
%
%INPUTS: c A scalar parameter > 0 that is the exponent of x in the PDF.
%        k A scalar parameter >0 that is the exponent of the denominator in
%          the PDF-1.
%   lambda A scalar parameter >0 that x is divided by in the PDF. This is
%          an extension from the original formulation of the distribution
%          in that x can be scaled. If this parameter is omitted or an
%          empty matrix is passed, then the default value of 1 is used.
%
%OUTPUTS: val The variance of the Burke type XII distribution.
%
%An expression for arbitrary noncentral moments of the distribution is
%given in [1]. This function just uses the identity that
%variance=E{x^2}-E{x}^2
%
%EXAMPLE:
%We verify the computed variance by comparing it to the sample variance.
% c=1.75;
% k=3;
% lambda=7;
% meanVal=BurrXIID.var(c,k,lambda)
% numSamp=1e6;
% sampMeanVal=var(BurrXIID.rand([numSamp,1],c,k,lambda))
%One will see that both values are about 9.2559.
%
%REFERENCES:
%[1] R. N. Rodriguez, "A guide to the Burr type XII distributions,"
%    Biometrica, vol. 64, no. 1, pp. 129-134, Apr. 1977.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<3||isempty(lambda))
       lambda=1; 
    end

    val=BurrXIID.kthMoment(2,c,k,lambda)-BurrXIID.kthMoment(1,c,k,lambda)^2;
end
    
function val=PDF(x,c,k,lambda)
%%PDF Evaluate the probability density function (PDF) of the Burr type XII
%     distribution.
%
%INPUTS: x The point or a matrix of points at which the PDF of the Burr
%          type XII distribution should be evaluated. The PDF is only
%          nonzero for positive x.
%        c A scalar parameter > 0 that is the exponent of x in the PDF.
%        k A scalar parameter >0 that is the exponent of the denominator in
%          the PDF-1.
%   lambda A scalar parameter >0 that x is divided by in the PDF. This is
%          an extension from the original formulation of the distribution
%          in that x can be scaled. If this parameter is omitted or an
%          empty matrix is passed, then the default value of 1 is used.
%
%OUTPUTS: val The value(s) of the PDF evaluated at the point(s) in x.
%
%The PDF is given in [1]. The incorporation of lambda is a simple change of
%variables.
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% c=1.75;
% k=3;
% lambda=7;
% numSamples=1e6;
% figure(1)
% clf
% histogram(BurrXIID.rand([numSamples,1],c,k,lambda),'Normalization','pdf')
% hold on
% numPoints=1000;
% x=linspace(0,20,numPoints);
% vals=BurrXIID.PDF(x,c,k,lambda);
% plot(x,vals,'linewidth',2)
% axis([0,20,0,0.2])
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%REFERENCES:
%[1] R. N. Rodriguez, "A guide to the Burr type XII distributions,"
%    Biometrica, vol. 64, no. 1, pp. 129-134, Apr. 1977.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.   

    if(nargin<4||isempty(lambda))
       lambda=1; 
    end

    val=(c*k/lambda)*(x/lambda).^(c-1)./(1+(x/lambda).^c).^(k+1);
    val(x<0)=0;
end
    
function prob=CDF(x,c,k,lambda)
%%CDF Evaluate the cumulative distribution function (CDF) of the Burr type
%     XII distribution.   
%
%INPUTS: x The point or a matrix of points at which the CDF of the Burr
%          type XII distribution should be evaluated. The CDF is only
%          nonzero for positive x.
%        c A scalar parameter > 0 that is the exponent of x in the PDF.
%        k A scalar parameter >0 that is the exponent of the denominator in
%          the PDF-1.
%   lambda A scalar parameter >0 that x is divided by in the PDF. This is
%          an extension from the original formulation of the distribution
%          in that x can be scaled. If this parameter is omitted or an
%          empty matrix is passed, then the default value of 1 is used.
%
%OUTPUTS: prob The value(s) of the CDF of the Burr type XII distribution
%              evaluated at x.
%
%The CDF is given in [1]. The incorporation of lambda is a simple change of
%variables.
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% c=1.75;
% k=3;
% lambda=7;
% x=5;
% numSamples=1e6;
% prob=BurrXIID.CDF(x,c,k,lambda)
% probSamp=mean(BurrXIID.rand([numSamples,1],c,k,lambda)<=x)
%One will see that both values are about 0.734.
%
%REFERENCES:
%[1] R. N. Rodriguez, "A guide to the Burr type XII distributions,"
%    Biometrica, vol. 64, no. 1, pp. 129-134, Apr. 1977.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.   

    if(nargin<4||isempty(lambda))
       lambda=1; 
    end

    prob=1-1./(1+(x/lambda).^c).^k;
    prob(x<=0)=0;
end
    
function x=invCDF(prob,c,k,lambda)
%%INVCDF Evaluate the inverse of the cumulative distribution function (CDF)
%        of the Burr Type XII distribution.
%
%INPUTS: prob The probability or probabilities (0<=prob<=1) at which the 
%             argument of the CDF is desired.
%           c A scalar parameter > 0 that is the exponent of x in the PDF.
%           k A scalar parameter >0 that is the exponent of the denominator
%             in the PDF-1.
%      lambda A scalar parameter >0 that x is divided by in the PDF. This
%             is an extension from the original formulation of the
%             distribution in that x can be scaled. If this parameter is 
%             omitted or an empty matrix is passed, then the default value
%             of 1 is used.
%
%OUTPUTS: x The argument(s) of the CDF that would give the probability or
%           probabilities in prob.
%   
%The CDF is given in [1]. The incorporation of lambda is a simple change of
%variables. Inversion of the CDF is straightforward.
%
%EXAMPLE:
%Here, we validate the inverse CDF by showing it to be the inverse of the
%CDF.
% c=1.75;
% k=3;
% lambda=7;
% x=5;
% xBack=BurrXIID.invCDF(BurrXIID.CDF(x,c,k,lambda),c,k,lambda)
%One will see that xBack is the same as x.
%
%REFERENCES:
%[1] R. N. Rodriguez, "A guide to the Burr type XII distributions,"
%    Biometrica, vol. 64, no. 1, pp. 129-134, Apr. 1977.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 
    
    if(nargin<4||isempty(lambda))
       lambda=1; 
    end

    x=lambda*nthroot(nthroot(1./(1-prob),k)-1,c);
end
    
function x=rand(N,c,k,lambda)
%%RAND Generate Burr type XII distributed random variables.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%        c A scalar parameter > 0 that is the exponent of x in the PDF.
%        k A scalar parameter >0 that is the exponent of the denominator in
%          the PDF-1.
%   lambda A scalar parameter >0 that x is divided by in the PDF. This is
%          an extension from the original formulation of the distribution
%          in that x can be scaled. If this parameter is omitted or an
%          empty matrix is passed, then the default value of 1 is used.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Burr type XII random variables.
%
%This function implements the inverse transformation method of Chapter 5.1
%of [1].
%
%REFERENCES:
%[1] S. M. Ross, Simulation, Ed. 4, Amsterdam: Elsevier, 2006.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C. 

    if(isscalar(N))
        dims=[N N];
    else
        dims=N;
    end
    
    if(nargin<4||isempty(lambda))
    	lambda=1; 
    end
    
    U=rand(dims);
    x=BurrXIID.invCDF(U,c,k,lambda);
end

function val=kthMoment(r,c,k,lambda)
%%KTHMOMENT Evaluate the rth noncentral moment of the Burr type XII
%           distribution.
%
%INPUTS: r The positive integer order of the moment. r<c*k for the moment
%          to exist. An error is thrown is the moment does not exist.
%        c A scalar parameter > 0 that is the exponent of x in the PDF.
%        k A scalar parameter >0 that is the exponent of the denominator in
%          the PDF-1.
%   lambda A scalar parameter >0 that x is divided by in the PDF. This is
%          an extension from the original formulation of the distribution
%          in that x can be scaled. If this parameter is omitted or an
%          empty matrix is passed, then the default value of 1 is used.
%
%OUTPUTS: The value fo the rth moment.
%
%The expression for arbitrary noncentral moments is given in terms of the
%beta function in [1]. A change of variables deals with lambda.
%
%REFERENCES:
%[1] R. N. Rodriguez, "A guide to the Burr type XII distributions,"
%    Biometrica, vol. 64, no. 1, pp. 129-134, Apr. 1977.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.  
    
    if(nargin<4||isempty(lambda))
       lambda=1; 
    end
    
    if(r>=c*k)
       error('The selected moment does not exist.'); 
    end

    val=k*(lambda)^r*beta(k-r/c,r/c+1);
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
