classdef WrappedNormalD
%%WRAPPEDNORMALD Functions to handle the scalar wrapped normal
%          distribution. The distribution is wrapped to a region of 2*pi in
%          length. The distribution is given in Chapter 2.2.6 of [1] and 
%Implemented methods are: circVar, PDF, CDF, trigMoment, params4TrigMoment,
%                         rand, wrappedNormalConvDist
%
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=circVar(sigma2)
%%CIRCVAR Obtain the circular variance of a particular wrapped normal
%         distribution.
%
%   sigma2 The variance of the linear normal distribution that is
%          wrapped.
%
%The circular variance is defined as 1-abs(r), where r is the mean
%resultant value of the second moment. The wrapped normal distribution is
%given in Chapter 2.2.6 of [1] and Proposition 2.1 in [1] discusses how to
%compute trigonometric moments of wrapped distributions.
%
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    val=1-exp(-sigma2/2);
end
    
function vals=PDF(x,mu,sigma2,maxNumTerms)
%%PDF Evaluate the PDF of the wrapped normal distribution.
%
%INPUTS: x A vector or matrix of points at which one wishes to evaluate the
%          wrapped normal distribution.
%       mu The mean of the linear normal distribution that is wrapped.
%   sigma2 The variance of the linear normal distribution that is
%          wrapped.
% maxNumTerms The maximum number of wrapped terms left and right of the
%          center point that should be used for the evaluation. If this
%          paramter is omitted or an empty matrix is passed, then the
%          default of 100 is used.
%
%OUTPUTS: vals The values of the wrapped normal PDF evaluated at the points
%              in x. This has the same dimensionality as x.
%
%The wrapped normal distribution is given in Chapter 2.2.6 of [1]. The PDF
%is Equation 2.2.14.
%
%EXAMPLE:
%Here, we validate the PDF by generating random samples and comparing the
%PDF plot with a histogram of the random samples.
% mu=4;
% sigma2=2;
% numRuns=1e6;
% numPoints=1000;
% 
% vals=wrappedNormalD.rand([numRuns,1],mu,sigma2);
% 
% figure(1)
% clf
% hold on
% histogram(vals,'normalization','pdf')
% x=linspace(-pi,pi,numPoints);
% PDFVals=WrappedNormalD.PDF(x,mu,sigma2);
% plot(x,PDFVals,'linewidth',2)
% h1=xlabel('x');
% h2=ylabel('PDF(x)');
% set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
% set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
%One will see that the histogram matches well with the plot.
%
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(nargin<4||isempty(maxNumTerms))
        maxNumTerms=100; 
    end
    
    N=numel(x);
    vals=zeros(size(x));
    
    expCoeff=-1/(2*sigma2);
    
    x=wrapRange(x-mu,0,2*pi);
    
    for curX=1:N
        xCur=x(curX);
        
        val=exp(xCur^2*expCoeff);
    
        xUpper=xCur;
        xLower=xCur;
        for curTerm=1:maxNumTerms
            xUpper=xUpper+2*pi;
            xLower=xLower-2*pi;
        
            oldVal=val;
            val=val+exp(xUpper^2*expCoeff)+exp(xLower^2*expCoeff);
        
            if(oldVal==val)
                break; 
            end
        end
        vals(curX)=val;
    end
    
    vals=vals/sqrt(2*pi*sigma2);
end

function vals=CDF(x,mu,sigma2,startingPoint,maxNumTerms)
%%CDF Evaluate the cumulative distribution function of the wrapped normal
%     distribution.
%
%INPUTS: x A vector or matrix of points at which one wishes to evaluate the
%          wrapped normal distribution.
%       mu The mean of the linear normal distribution that is wrapped.
%   sigma2 The variance of the linear normal distribution that is
%          wrapped.
% startingPoint This is the point in the region from which integration is
%          started. Note that x=startingPoint+2*pi is equivalent to
%          x=startingPoint and will return a CDF value of 0 (not 1).
% maxNumTerms The maximum number of wrapped terms left and right of the
%          center point that should be used for the evaluation. If this
%          paramter is omitted or an empty matrix is passed, then the
%          default of 100 is used.
%
%OUTPUTS: vals The values of the wrapped normal CDF evaluated at the points
%              in x. This has the same dimensionality as x.
%
%The CDF comes from integration the wrapped normal PDF< which is given in
%Chapter 2.2.6 of [1]. Note that the CDF of a circular distribution is
%periodic.
%
%EXAMPLE:
%We validate the CDF value by comparing it to a value computed from random
%samples.
% mu=4;
% sigma2=2;
% x=2;
% numSamples=1e6;
% prob=WrappedNormalD.CDF(x,mu,sigma2)
% probSamp=mean(WrappedNormalD.rand([numSamples,1],mu,sigma2)<=x)
%One will see that both values are about 0.8056.
%   
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    if(nargin<5||isempty(maxNumTerms))
        maxNumTerms=100; 
    end
    
    if(nargin<4||isempty(startingPoint))
       startingPoint=-pi; 
    end

    N=numel(x);
    vals=zeros(size(x));
    
    startingPoint=mod(startingPoint,2*pi);
    mu=mod(mu,2*pi);
    for curX=1:N
        xCur=mod(x(curX),2*pi);

        val=GaussianCDFOverRegion(startingPoint,xCur,mu,sigma2);
        for k=1:maxNumTerms
            valOld=val;

            val=val+GaussianCDFOverRegion(startingPoint+2*pi*k,xCur+2*pi*k,mu,sigma2)+GaussianCDFOverRegion(startingPoint-2*pi*k,xCur-2*pi*k,mu,sigma2);
        
            if(val==valOld)
                break;
            end
        end

        if(xCur<startingPoint)
            val=1+val;
        end
        vals(curX)=val;
    end
end

function [rho,theta,R]=trigMoment(n,mu,sigma2)
%%TRIGMOMENT Compute the nth trigonometric moment of the wrapped normal
%            distribution. A direction theta on the unit circle can be
%            modeled as a complex quantity having unit magnitude as
%            exp(1j*theta). The nth raw complex trigonometric moment is the
%            expected value of exp(1j*n*theta), where the integral for the
%            expected value is taken over any interval of length 2*pi.
%
%INPUTS: n The order of the moment desired. This is >=1.
%       mu The mean of the linear normal distribution that is wrapped.
%   sigma2 The variance of the linear normal distribution that is
%          wrapped.
%
%OUTPUTS: rho The (complex) mean resultant value for the nth moment. Note
%             that rho=R*exp(1j*theta).
%       theta The (real) trigonometric mean angle in radians for the nth
%             moment. This is between -pi and pi.
%           R The (real) mean resultant length for the nth moment.
%
%Proposition 2.1 in [1] discusses how to compute trigonometric moments of
%wrapped distributions.
%
%EXAMPLE:
%We validate the third trigonometric moment value by comparing it to a
%value computed from random samples.
% mu=4;
% sigma2=0.8;
% n=3;
% numSamples=1e7;
% vals=WrappedNormalD.rand([numSamples,1],mu,sigma2);
% rhoSamp=mean(exp(1j*vals*n))
% rho=WrappedNormalD.trigMoment(n,mu,sigma2)
%One will see that both values are about 0.0231-0.0147i.
%
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C

rho=exp(1i*n*mu-n^2*sigma2/2);
if(nargout>1)
    theta=angle(rho);
    R=abs(rho);
end
end

function [mu,sigma2]=params4TrigMoment(rho)
%%PARAMS4TRIGMOMENT Obtain the parameters (mu and Sigma2) of the wrapped
%           normal distribution such that the wrapped normal distribution
%           matches the given complex mean resultant length for the first
%           moment.
%
%INPUTS: rho The (complex) mean resultant value for the first moment. Note
%            that abs(rho)<=1.
%
%OUTPUTS: mu, sigma2 The mean and variance of the linear normal
%                    distribution that is wrapped. mu is in the range -pi
%                    to pi.
%
%Proposition 2.1 in [1] discusses how to compute trigonometric moments of
%wrapped distributions. We just solve the inverse problem here.
%
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C

mu=atan2(imag(rho),real(rho));
sigma2=-2*log(abs(rho)); 
end

function vals=rand(N,mu,sigma2)
%%RAND Generate wrapped normal random variables.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%       mu The mean of the linear normal distribution that is wrapped.
%   sigma2 The variance of the linear normal distribution that is
%          wrapped.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated wrapped normal random variables. The values are
%              wrapped to the region -pi to pi.
%
%This function just generated normal random variables and wraps then to the
%region from -pi to pi usign wrapRange.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C    

    vals=wrapRange(mu+sqrt(sigma2)*randn(N),-pi,pi);
end


function [mu,sigma2]=wrappedNormalConvDist(mu1,sigma21,mu2,sigma22)
%WRAPPEDNORMALCONVDIST Compute the distribution from convolving two wrapped
%         normal random distributions. As noted in Chapter 2.2.6 of [1],
%         the result is also a wrapped normal distribution. This function
%         returns the parameters of that distribution.
%
%INPUTS: mu1, sigma21 The mean and variance of the linear normal
%                      distribution underlying the first wrapped normal
%                      distribution.
%        mu1, sigma21 The mean and variance of the linear normal
%                      distribution underlying the second wrapped normal
%                      distribution.
%
%OUTPUTS: mu, sigma2 The mean and variance of the wrapped normal
%                    distribution resulting from convolving the two input
%                    wrapped normal distributions. mu us wrapped to the
%                    range from -pi to pi.
%   
%REFERENCES:
%[1] S. R. Jammalamadaka and A. SenGupta, Topics in Circular Statistics.
%    Singapore: World Scientific, 2001.
%
%April 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.   
    
    mu=wrapRange(mu1+mu2,-pi,pi);
    sigma2=sigma21+sigma22;
end

end
end

function val=GaussianCDFOverRegion(startPos,endPos,mu,sigma2)
    denom=sqrt(2*sigma2);
    val=(1/2)*(erf((mu-startPos)/denom)-erf((mu-endPos)/denom));
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
