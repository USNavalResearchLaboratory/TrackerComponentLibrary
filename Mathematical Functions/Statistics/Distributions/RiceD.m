classdef RiceD
%%RICED Functions to handle the Rice distribution.
%Implemented methods are: mean, var, PDF, CDF, rand
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    
function val=mean(s,sigma)
%%MEAN Obtain the mean of the Rice distribution for given noncentrality
%      and scale parameters.
%
%INPUTS: s The noncentrality parameter of the distribution.
%    sigma The scale parameter of the distribution.
%
%OUTPUTS: val The mean of the Rice distribution.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    param=-s^2/(2*sigma^2);
    val=sigma*sqrt(pi/2)*exp(param/2)*((1-param)*besseli(0,-param/2)-param*besseli(1,-param/2));
end

function val=var(s,sigma)
%%VAR Obtain the variance of the Rice distribution for given noncentrality
%     and scale parameters.
%
%INPUTS: s The noncentrality parameter of the distribution.
%    sigma The scale parameter of the distribution.
%
%OUTPUTS: val The variance of the Rice distribution.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
    param=-s^2/(2*sigma^2);
    L=exp(param/2)*((1-param)*besseli(0,-param/2)-param*besseli(1,-param/2));
    val=2*sigma^2+s^2-pi*(sigma^2/2)*L^2;
end 
    
function val=PDF(x,s,sigma)
%%PDF Evaluate the Rice probability distribution function at one or more
%     desired points.
%
%INPUTS: x The point(s) at which the Rice PDF is to be evaluated. Note that
%          x>=0.
%        s The noncentrality parameter of the distribution.
%    sigma The scale parameter of the distribution.
%
%OUTPUTS: val The value(s) of the Rice PDF with parameters s and sigma
%             evaluated at x.
%
%The PDF of the Rice distribution is given in Chapter 2.1.4 of [1].
%
%REFERENCES:
%[1] J. G. Proakis, Digital Communications. Ed. 4, Boston, MA:
%    McGraw Hill, 2001.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    %We are evaluating 
    %val=(x/sigma^2).*exp(-(x.^2+s^2)/(2*sigma^2)).*besseli(0,x*s/sigma^2);
    %However, if x*s/sigma^2 is large, then this function could return
    %NaNs. large doesn't have to be all that large. For example, s=30,
    %sigma=1, x=100. means that ratio is 3000. Thus, we use besseli with
    %the third argument to get a scaled value, then take the logarithm, do
    %everything else in the logarithmic domain and take the exponent to
    %undo the logarithm.
    
    besselArg=x*s/sigma^2;
    val=exp(log(besseli(0,besselArg,1))+besselArg-(x.^2+s^2)/(2*sigma^2)+log(x/sigma^2));
end

function val=CDF(x,s,sigma)
%%CDF Evaluate the cumulative distribution function of the Rice
%     distribution at desired points.
%
%INPUTS: x The point(s) at which the Rice CDF is to be evaluated. Note that
%          x>=0.
%        s The noncentrality parameter of the distribution.
%    sigma The scale parameter of the distribution.
%
%OUTPUTS: prob The value(s) of the CDF of the Rice distribution with
%              parameters k and theta evaluated at x.
%
%The CDF of the Rice distribution can be expressed in terms of Marcum's Q
%function as shown in Chapter 2.1.4 of [1].
%
%REFERENCES:
%[1] J. G. Proakis, Digital Communications. Ed. 4, Boston, MA:
%    McGraw Hill, 2001.
%
%February 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
    numPoints=length(x(:));
    val=zeros(size(x));
    
    for curPoint=1:numPoints
        val(curPoint)=1-MarcumQ(1,s/sigma,x(curPoint)/sigma);
    end
end


function vals=rand(N,s,sigma)
%%RAND Generate Rice-distributed random variables with the given
%      parameters.
%
%INPUTS: N If N is a scalar, then rand returns an NXN matrix of random
%          variables. If N=[M,N1] is a two-element row vector, then rand
%          returns an MXN1 matrix of random variables.
%        s The noncentrality parameter of the distribution.
%    sigma The scale parameter of the distribution.
%
%OUTPUTS: vals A matrix whose dimensions are determined by N of the
%              generated Rice random variables.
%
%This generates Rice distributed random variables by transforming normally
%distributed random variables using the identity given in Chapter 2.1.4 of
%[1].
%
%REFERENCES:
%[1] J. G. Proakis, Digital Communications. Ed. 4, Boston, MA:
%    McGraw Hill, 2001.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

    if(isscalar(N))
        dims=[N, N];
    else
        dims=N;
    end

    m=s/sqrt(2);

    X=sigma*randn(dims)+m;
    Y=sigma*randn(dims)+m;

    vals=sqrt(X.^2+Y.^2);
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
