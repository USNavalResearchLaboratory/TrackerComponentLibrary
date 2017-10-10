function f=ChebyshevPolySynth(tau,a,tauStart,tauEnd)
%%CHEBYSHEVPOLYSYNTH Evaluate a weighted sum of Chebyshev polynomials of
%                    the first kind given the weighting coefficients. This
%                    is generally done to perform interpolation. This
%                    function is more efficient than using the function
%                    ChebyshevPoly to evaluate each of the Chebyshev
%                    polynomials individually.
%
%INPUTS: tau An NX1 or 1XN vector of values from tauStart to tauEnd, or
%            from -1 to 1 if tauStart and tauEnd are omitted, where one
%            wishes to evaluate the weighted Chebyshev polynomial series.
%          a An (n+1)X1 or 1X(n+1) vector of the coefficients for the
%            Chebyshev polynomials from degree 0 to degree n.
% tauStart,tauEnd The possible range of the inputs. If omitted, a  range of
%            -1 to 1 is assumed --the normal range for Chebyshev
%            polynomials. The option for mapping to a wider range is useful
%            when using Chebyshev polynomials for interpolation. Note that
%            the such polynomials are generally not useful for
%            interpolating much outside of the valid range.
%
%OUTPUTS: f An NX1 vector of the weighted Chebyshev polynomial series sum
%           evaluated at all of the points in tau.
%
%The algorithm used is the Clenshaw method, which is described in Chapter
%3.3.3 of [1].
%
%REFERENCES:
%[1] O. Montenbruck and E. Gill, Satellite Orbits, 4th ed. Heidelberg:
%    Springer, 2012.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(a)-1;
%Make tau a column vector.
tau=tau(:);

%Map the tau values to the -1 to 1 range, if necessary.
if(nargin>2)
    tau=(tau-0.5*(tauStart+tauEnd))/(0.5*(tauEnd-tauStart));
end

fnp1=0;
fnp2=0;

for curN=n:-1:1
    f=2*tau.*fnp1-fnp2+a(curN+1);
    fnp2=fnp1;
    fnp1=f;
end

f=tau.*fnp1-fnp2+a(1);
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
