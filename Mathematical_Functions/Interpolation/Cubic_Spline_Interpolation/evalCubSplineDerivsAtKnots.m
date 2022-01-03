function derivVals=evalCubSplineDerivsAtKnots(C,tau,derivNum)
%%EVALCUBSPLINEDERIVSATKNOTS Given the output of the fitCubSpline function,
%   obtain the first, second or third derivatives evaluated at the knots.
%   This is faster than calling evalCubSpline for all of the knots.
%
%INPUTS: C The 4XN matrix of coefficients returned by the fitCubSpline
%          function.
%      tau The 1XN or NX1 vector of independent variable points fitted in
%          the fitCubSpline function that produced C (breakpoints).
%  derivNum The number of derivatives desired. This is a value from 0-3.
%           The default if omitted or an empty matrix is passed is 1.
%
%OUTPUTS: derivVals A 1:N vector of the derivNum derivative (N is the
%                   number of knots used to create pp).
%
%The derivatives for points 1:(N-1) are given directly in the coefficients
%in C. The final point must be interpolated from the coefficients. The
%interpolation is described in Chapter 4 of [1]. More details on the
%formulae are in them comments to the function evalCubSpline.
%
%REFERENCES:
%[1] C. de Boor, A Practical Guide to Splines. New York: Springer-Verlag,
%    1978.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(derivNum))
    derivNum=1;
end

N=length(tau);

switch(derivNum)
    case 0
        derivVals=C(1,:);
    case 1
        diff=tau(N)-tau(N-1);
        endDeriv=C(2,N-1)+C(3,N-1)*diff+C(4,N-1)*diff^2/2;
        derivVals=[C(2,1:(N-1)),endDeriv];
    case 2
        diff=tau(N)-tau(N-1);
        endDeriv=C(3,N-1)+C(4,N-1)*diff;
        derivVals=[C(3,1:(N-1)),endDeriv];
    case 3
        derivVals=[C(4,1:(N-1)),C(4,(N-1))];
    otherwise
        error('Invalid derivative number specified.')
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
