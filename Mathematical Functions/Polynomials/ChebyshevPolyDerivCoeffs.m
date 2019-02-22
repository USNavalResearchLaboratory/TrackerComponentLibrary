function b=ChebyshevPolyDerivCoeffs(a,tauStart,tauEnd)
%%CHEBYSHEVPOLYDERIVCOEFFS Given the coefficients to a series of Chebyshev
%                   polynomials of the first kind, usually used for
%                   interpolation, that one could evaluate using the
%                   ChebyshevPolySynth function, compute the coefficients
%                   for the derivative of the function.
%
%INPUTS: a A 1X(n+1) or (n+1)X1 vector of coefficients for a weighted sum
%          of Chebyshev polynomials of order 0 to n.
% tauStart,tauEnd The possible range of the inputs to the Chebyshev
%          functions. This is used if the functions are transformed to map
%          to a range larger than their standard -1 to 1. If omitted, a
%          range of -1 to 1 is assumed --the normal range for Chebyshev
%          polynomials. The option for mapping to a wider range is useful
%          when using Chebyshev polynomials for interpolation. Note that
%          such polynomials are generally not useful for interpolating far
%          outside of the valid range.
%
%INPUTS: b An nX1 vector of coefficients for a weighted sum of Chebyshev
%          polynomials of order 0 to n-1 that can be used in the
%          ChebyshevPolySynth function to find the derivative of the
%          function given by a. If a is 1X1, then b=0 will be returned.
%
%The derivative of a weighted Chebyshev polynomial series is another
%weighted Chebyshev polynomial series, as given in [1]. The derivative of
%the Chebyshev polynomial series is exact, as noted in the above, it comes
%from a simplification of term-by-term differentiation. However, it is not
%necessarily as accurate as one might hope when the polynomials are used
%for interpolation. Though the above notes a number of issues in
%interpolated derivatives, one would expect such problems to be lessened if
%the derivatives themselves were used when designing the interpolation
%routine. For example, in [2], the Chebyshev polynomials are fit to points
%and the first derivatives.
%
%The opposite of this function is ChebyshevPolyIntCoeffs, though it cannot
%be strictly used as an inverse, because it assumes the additive
%integration constant is zero.
%
%REFERENCES:
%[1] K. S. Breuer and R. M. Everson, "On the errors incurred calculating
%    derivatives using Chebyshev polynomials," Journal of Computational
%    Physics, vol. 99, no. 1, pp. 56-67, Mar. 1992.
%[2] X. X. Newhall, "Numerical representation of planetary ephemerides,"
%    Celestial Mechanics, vol. 45, no. 1-3, pp. 305-310, 1989.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(a)-1;

%Allocate space.
b=zeros(n+1,1);

%Set the default range.
if(nargin<2)
    tauEnd=1;
    tauStart=-1;
end

%Set the last elements
if(n==0)
    b=0;
elseif(n==1)
    k=1;
    b(k+1)=0;
    k=0;
    b(k+1)=a(k+1+1);
else
    k=n;
    b(k+1)=0;
    k=n-1;
    b(k+1)=2*(k+1)*a(k+1+1);
    %Perform recursion for the rest of the elements.
    for k=(n-2):-1:1
       b(k+1)=b(k+2+1)+2*(k+1)*a(k+1+1); 
    end

    k=0;
    b(k+1)=b(k+2+1)/2+a(k+1+1);
end

%Get rid of the last (zero) element, unless it would make b an empty
%matrix.
if(n>0)
b=b(1:end-1);
end

%The constant to handle a mapping from -1->1 to some other range of
%parameters.
b=b*(2/(tauEnd-tauStart));

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
