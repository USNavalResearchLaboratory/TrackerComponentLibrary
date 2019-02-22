function [dT,T]=ChebyshevPolyDeriv(tau,n,tauStart,tauEnd)
%%CHEBYSHEVPOLYDERIV Compute the derivatives of Chebyshev polynomials of
%                    the first kind from order 0 to order n evaluated at
%                    the points given in tau. The explicit derivatives can
%                    be useful for least-squares fitting of a weighted
%                    Chebyshev polynomial to a given function given the
%                    function evaluated at certain points along with the
%                    derivatives of the function at those points.
%
%INPUTS: tau An NX1 or 1XN vector of values from tauStart to tauEnd, or
%            from -1 to 1 if tauStart and tauEnd are omitted, where one
%            wishes to evaluate the derivatives of the Chebyshev
%            polynomials.
%          n The non-negative integer maximum order of the Chebyshev
%            polynomials evaluated.
% tauStart,tauEnd The possible range of the inputs. If omitted, a range of
%            -1 to 1 is assumed --the normal range for Chebyshev
%            polynomials. The option for mapping to a wider range is useful
%            when using Chebyshev polynomials for interpolation. Note that
%            such polynomials are generally not useful for interpolating
%            far outside of the valid range.
%
%OUTPUTS: dT An (n+1)XN matrix of the derivatives of the Chebyshev
%            polynomials from order 0 to n evaluated at each of the values
%            of tau. T(i,j) is the derivative of the (i-1)th order
%            Chebyshev polynomial evaluated at tau(j), taking into account
%            that the function has been mapped to the range
%            (tauStart,tauEnd) if necessary.
%          T Since the Chebyshev functions have to be computed to find the
%            derivatives, an (n+1)XN matrix of the Chebyshev functions can
%            also be returned, if desired.
%
%The recursion for the derivatives come from simply differentiating the
%recursion formula for the Chebyshev polynomials themselves, which is in
%[1] and in Chapter 3.3.3 of [2].
%
%REFERENCES:
%[1] X. X. Newhall, "Numerical representation of planetary ephemerides,"
%    Celestial Mechanics, vol. 45, no. 1-3, pp. 305-310, 1989.
%[2] O. Montenbruck and E. Gill, Satellite Orbits, 4th ed. Heidelberg:
%    Springer, 2012.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=length(tau);
tau=tau(:)';%Make tau into a column vector.
%Map the tau values to the -1 to 1 range, if necessary.
if(nargin>2)
    tau=(tau-0.5*(tauStart+tauEnd))/(0.5*(tauEnd-tauStart));
end
%Get the Chebyshev polynomial function values.
T=ChebyshevPoly(tau,n);
dT=zeros(n+1,N);

%The zeroth-order derivative is zero.
dT(0+1,:)=0;
if(n>0)
    dT(1+1,:)=1;
    for k=2:n
        dT(k+1,:)=2*tau.*dT(k-1+1,:)+2*T(k-1+1,:)-dT(k-2+1,:);
    end
end

%The constant to handle a mapping from -1->1 to some other range of
%parameters.
if(nargin>2)
    dT=dT*(2/(tauEnd-tauStart));
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
