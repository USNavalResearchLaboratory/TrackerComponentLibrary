function T=ChebyshevPoly1(tau,n,tauStart,tauEnd)
%%CHEBYSHEVPOLY1 Evaluate a Chebyshev polynomial of the first kind of a
%                single given order at the specified points. If one wishes
%                to evaluate all Chebyshev polynomials from order 0 to n
%                then use the function ChebyshevPoly. This function will
%                also evaluate a generalized Chebyshev polynomial in which
%                case tau is not mapped to the region [-1,1].
%
%INPUTS: tau An NX1 or 1XN vector of values from tauStart to tauEnd, or
%            from -1 to 1 if tauStart and tauEnd are omitted, where one
%            wishes to evaluate the Chebyshev polynomials.
%          n The non-negative integer order of the Chebyshev polynomial
%            evaluated.
% tauStart,tauEnd The possible range of the inputs. Unless the input is
%            unbounded as described below, given this range, the actual
%            input to the polynomials is scaled to
%            tau=(tau-0.5*(tauStart+tauEnd))/(0.5*(tauEnd-tauStart));
%            If omitted, a  range of -1 to 1 is assumed --the normal range
%            for Chebyshev polynomials. The option for mapping to a wider
%            range is useful when using Chebyshev polynomials for
%            interpolation. Note that such polynomials are generally not
%            useful for interpolating far outside of the valid range. If
%            one wishes to use generalized Chebyshev polynomials that are
%            also defined outside of the range of [-1,1], then use
%            tauStart=-Inf,tauEnd=Inf. The generalized polynomials are
%            described below.
%
%OUTPUTS: T An NX1 or 1XN vector of the Chebyshev polynomial of order n
%           evaluated at the given point(s).
%
%Where as the function ChebyshevPoly uses a recursion to evaluate the
%Chebyshev polynomials of order 0 to n at the points, that is not efficient
%if one wishes to evaluate a single Chebyshev polynomial when n is large.
%This function uses the identity from [1] that Chebyshev polynomials can be
%expressed in terms of trigonometric functions.
%
%If tauStart=-Inf and TauEnd=Inf, then tau is not scaled to the range -1 to
%1 based on the bounds. Instead, the extension of Chebyshev polynomials to
%all real numbers as given in Chapter 3.4.2 of [2] is used. Outside of the
%range of -1 to 1, cosine terms are replaced with cosh terms in the
%evaluation and there is a (-1)^n coefficient out front for tau<-1.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Chebyshev Polynomial of the First Kind." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html
%[2] H. L. Van Trees, Optimum Array Processing. New York: Wiley-
%    Interscience, 2002.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Map the tau values to the -1 to 1 range, if necessary.
if(nargin>2)
    if(tauStart==-Inf&&tauEnd==Inf)
        %In this case, we are evaluating a generalized Chebyshev
        %polynomial.
        sel1=abs(tau)<=1;
        sel2=tau>1;
        sel3=tau<-1;
        
        T=zeros(size(tau));
        T(sel1)=cos(n*acos(tau(sel1)));
        T(sel2)=cosh(n*acosh(tau(sel2)));
        T(sel3)=(-1)^n*cosh(n*acosh(tau(sel3)));
        return;
    else
        tau=(tau-0.5*(tauStart+tauEnd))/(0.5*(tauEnd-tauStart));
    end
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
