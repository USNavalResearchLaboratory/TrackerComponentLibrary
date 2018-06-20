function a=ChebyshevPolyFit(f,n,tauStart,tauEnd)
%%CHEBYSHEVPOLYFIT Get weights to fit a series of Chebyshev polynomials of
%                  the first kind to a function f over the range of values
%                  tauStart,tauEnd. The coefficients can be used with the
%                  function ChebyshevPolySynth to interpolate the
%                  function. Fitted Chebyshev polynomials are a good
%                  approximation to a minimax polynomial fit. In such an
%                  instance, one generally chooses n notably larger than
%                  desired, and then truncates the resulting coefficient
%                  set.
%
%INPUTS: f A function handle to a scalar function that can be passed
%          vectors of points at which it is to be evaluated for
%          interpolation. The function will not be evaluated at the
%          endpoints.
%        n The maximum order of the polynoimals to generate. There will be
%          n+1 returned coefficients, because there is a 0 order
%          polynomial.
% tauStart,tauEnd The range of the f over which fiting for interpolation is
%          to be taken. If omitted, a range of -1 to 1 is assumed. This is
%          the normal range for Chebyshev polynomials.
%
%OUTPUTS: a An (n+1)X1 or 1X(n+1) vector of the coefficients for the
%           Chebyshev polynomials from degree 0 to degree n. This can be
%           used with the function ChebyshevPolySynth for interpolation.
%
%The coefficients are deterministically found by evaluating the function f
%at a predictable set of points based on the Chebyshev Approximation
%Formula that is mentioned in [1]. At each of the points chosen for
%evaluation, which are the zeros of the highest-order Chebyshev polynomial,
%the function is matched exactly.
%
%The discrete costine transform was used to implement the sum in a
%computationally efficient manner.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Chebyshev Approximation Formula." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/ChebyshevApproximationFormula.html
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(tauStart))
    tauStart=-1;
    tauEnd=1;
end

%The number of points is one more than the maximum order.
N=n+1;

k=1:N;
cosVals=cos(pi*(k-0.5)/N);
%The cosine values must be mapped to the range tauStart,tauEnd.
x=0.5*(tauStart+tauEnd+cosVals*(tauEnd-tauStart));
fVals=f(x);

a=(1/N)*discSinCosTrans(fVals,'CIIe');
a(1)=0.5*a(1);

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
