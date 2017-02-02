function T=ChebyshevPoly1(tau,n,tauStart,tauEnd)
%%CHEBYSHEVPOLY1 Evaluate a Chebyshev polynomial of the first kind of a
%                single given order at the specified points. If one wishes
%                to evaluate all Chebyshev polynomials from order 0 to n
%                then use the function ChebyshevPoly.
%
%INPUTS: tau     An NX1 or 1XN vector of values from tauStart to tauEnd, or
%                from -1 to 1 if tauStart and tauEnd are omitted, where one
%                wishes to evaluate the Chebyshev polynomials.
%        n       The non-negative integer order of the Chebyshev polynomial
%                evaluated.
%
%tauStart,tauEnd The possible range of the inputs. If omitted, a  range of
%               -1 to 1 is assumed --the normal range for Chebyshev
%                polynomials. The option for mapping to a wider range is
%                useful when using Chebyshev polynomials for interpolation.
%                Note that the such polynomials are generally not useful
%                for interpolating much outside of the valid range.
%
%OUTPUTS: T      An NX1 or 1XN vector of the Chebyshev polynomial of order
%                n evaluated at the given point(s).
%
%Where as the function ChebyshevPoly uses a recursion to evaluate the
%Chebyshev polynomials of order 0 to n at the points, that is not efficient
%if one wishes to evaluate a single Chebyshev polynomial when n is large.
%This function uses the identity from [1] that Chebyshev polynomials can be
%expressed in terms of trigonometric functions.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Chebyshev Polynomial of the First Kind." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Map the tau values to the -1 to 1 range, if necessary.
if(nargin>2)
    tau=(tau-0.5*(tauStart+tauEnd))/(0.5*(tauEnd-tauStart));
end

T=cos(n*acos(tau));

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
