function [intEst,absErrEst,exitCode]=integral1DAdaptiveCC(func,bounds,RelTol,AbsTol,NPowMin,NPowMax)
%%INTEGRAL1DADAPTIVECC Perform adaptive numerical integration over a
%           function that takes a scalar input and returns a scalar or
%           vector quantity. The integration is performed using
%           Clenshaw-Curtis points. The function ClenshawCurtisPoints1D is
%           used to get points of two orders and the difference in the
%           integral values from each is the error estimate. If the error
%           is too high, then the next higher pair of orders is used (the
%           next one being double the previous order). Note
%           that Clenshaw-Curtis points of length 2*N have all the
%           quadrature points of length N embedded, so by using powers of
%           2 for the lengths, we can reuse previously computed function
%           values.
%
%INPUTS: f The handle to the function. f(x) must be able to take a 1XN
%          vector of points and return the numDimXN set of function values
%          at those points.
%   bounds A 2X1 or 1X2 vector such that bounds(1) is the lower bound of
%          integration and bounds(2) is the upper bound of integration. it
%          is not required that bounds(2)>bounds(1). Integration bounds
%          must be finite. The default if omitted or an empty matrix is
%          passed is [-1;1].
%   RelTol The maximum relative error tolerance allowed. If omitted or an
%          empty matrix is passed, the default value of 1e-8 is used. When
%          f returns a vector, this is applied to each element separately.
%   AbsTol The absolute error tolerance allowed. If omitted or an empty
%          matrix is passed, the default value of 1e-11 is used. When
%          f returns a vector, this is applied to each element separately.
%  NPowMin The minimum order of the Chebyshev polynomials to use is
%          2^NPowMin. NPowMin is an integer>=1. The default if omitted or
%          an empty matrix is passed is 2.
%  NPowMax The maximum order of the Chebyshev polynomials that will be
%          considered is 2^NPowMax. The default if omitted or
%          an empty matrix is passed is 11, referring to order 2048.
%          NPowMax must be greater than NPowMin.
%
%OUTPUTS: intEst The estimate of the integral.
%      absErrEst An estimate of the error of the integral. This is the
%                difference of the values at two different orders.
%       exitCode A value indicating the status of the algorithm on
%                termination. Possible values are
%                0 Termination occurred due to the absolute or relative
%                  error tolerances being fulfilled.
%                1 Termination occurred due to reaching NPowMax.
%                2 A singularity of a NaN was encountered during a step;
%                  the returned values of intEst and totalError are from
%                  the previous step and do not meet the error criterion.
%                3 A singularity or NaN was encountered. The values of
%                  intEst and totalError are most likely invalid.
%
%The initial Chebyshev interpolation attempted will be order 2^NPowMin and
%an approximation of order 2^(NPowMin+1) will be used to get the error
%estimate. After that, the order of the pairs attempted is increased until
%convergence or the last error estimate was obtained having order NPowMax
%without convergence.
%
%EXAMPLE:
%Here, we compare the integral of a nonlinear function that is rather peaky
%in the beginning with the Chebyshev approximated integral.
% f=@(x)(x.*exp(-x.^2-x/2));
% bounds=[0;16];
% [intEst,absErrEst,exitCode]=integral1DAdaptiveCC(f,bounds);
% explicitSol=(1/8)*(4-4/exp(264)+exp(1/16)*sqrt(pi)*(erf(1/4)-erf(65/4)));
% relErr=abs(intEst-explicitSol)/abs(explicitSol)
%The relative error will be on the order of 1e-16, indicating a good
%agreement between the explicit an numerically determined integrals.
%
%May 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(NPowMax))
    NPowMax=11;
end

if(nargin<5||isempty(NPowMin))
    NPowMin=2;
end

if(NPowMax<=NPowMin)
    error('It is required that NPowMax>NPowMin.');
end

if(nargin<4||isempty(AbsTol))
    AbsTol=1e-11;
end

if(nargin<3||isempty(RelTol))
    RelTol=1e-8;
end

if(any(~isfinite(bounds(:))))
    error('The bounds must be finite. Perform a change of variables to handle infinite bounds.');
end

if(nargin<2||isempty(bounds)||(bounds(1)==-1&&bounds(2)==1))
    %Integration bounds are -1 to 1.
    f=func;
    df=1;
else
    %Integration bounds are given in bounds. We have to map them to (-1,1).
    LB=bounds(1);
    UB=bounds(2);
    f=@(x)(func(((UB-LB)/2).*x+(UB+LB)/2));
    df=(UB-LB)/2;
end
NCur=NPowMin;
[xi,w]=ClenshawCurtisPoints1D(2^NCur);
fLo=f(xi);
numDim=size(fLo,1);
%The integral value with the lower polynomial order approximation.
ILow=sum(fLo.*w',2);

%If a non-finite term was encountered.
if(~isfinite(ILow))
    absErrEst=Inf;
    exitCode=3;
    intEst=ILow*df;
    return;
end

absErrEstPrev=Inf;
NCur=NCur+1;
while(NCur<=NPowMax)
    highPow=2^NCur;
    [xi,w]=ClenshawCurtisPoints1D(2^NCur);
    %The Clenshaw-Curtis points are nested, so we do not need to reevaluate
    %the function at the points it was already evaluated at (every other
    %point).
    numPts=highPow+1;
    fHigh=zeros(numDim,numPts);
    fHigh(:,1+(0:2:highPow))=fLo;
    fHigh(:,1+(1:2:(highPow-1)))=f(xi(1+(1:2:(highPow-1))));

    IHigh=sum(fHigh.*w',2);

    %If a non-finite term was encountered.
    if(~isfinite(ILow))
        absErrEst=absErrEstPrev;
        exitCode=2;
        intEst=ILow*df;
        return;
    end
    
    %Check convergence
    absErrEst=abs(IHigh-ILow);
    if(all(absErrEst<=AbsTol)||all(absErrEst<=RelTol*abs(IHigh)))
        exitCode=0;
        intEst=IHigh*df;
        return
    end

    absErrEstPrev=absErrEst;
    ILow=IHigh;
    fLo=fHigh;
    NCur=NCur+1;
end

%If the maximum power was reached.
exitCode=1;
intEst=IHigh*df;

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
