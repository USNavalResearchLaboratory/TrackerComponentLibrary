function b=ChebyshevPolyIntCoeffs(a,tauStart,tauEnd)
%%CHEBYSHEVPOLYINTCOEFFS Given the coefficients to a series of Chebyshev
%                   polynomials of the first kind, usually used for
%                   interpolation, that one could evaluate using the
%                   ChebyshevPolySynth function, compute the coefficients
%                   for the indefinite integral of the function.
%
%INPUTS: a A 1X(n+1) or (n+1)X1 vector of coefficients for a weighted sum
%          of Chebyshev polynomials of order 0 to n.
%tauStart,tauEnd The possible range of the inputs to the Chebyshev
%          unctions. This is used if the functions are transformed to map
%          to a range larger than their standard -1 to 1. If omitted, a
%          range of -1 to 1 is assumed --the normal range for Chebyshev
%          polynomials. The option for mapping to a wider range is useful
%          when using Chebyshev polynomials for interpolation. Note that
%          such polynomials are generally not useful for interpolating far
%          outside of the valid range.
%
%INPUTS: b An (n+2)X1 vector of coefficients for a weighted sum of
%          Chebyshev polynomials of order 0 to n+1 that can be used in the
%          ChebyshevPolySynth function to find the indefinte integral of
%          the function given by a (with 0 additive constant).
%
%The derivative of a weighted Chebyshev polynomial series is another
%weighted Chebyshev polynomial series, as given in [1] and is derived by
%element-by-element differentiation. The integral can be derived in a
%similar manner.
%
%The opposite of this function is ChebyshevPolyDerivCoeffs.
%
%REFERENCES:
%[1] K. S. Breuer and R. M. Everson, "On the errors incurred calculating
%    derivatives using Chebyshev polynomials," Journal of Computational
%    Physics, vol. 99, no. 1, pp. 56-67, Mar. 1992.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The order of a
n=length(a)-1;

%Allocate space.
b=zeros(n+2,1);

%Set the default range.
if(nargin<2)
    tauEnd=1;
    tauStart=-1;
end
 
%The constant to normalize the coefficients to a given range
normConst=(tauEnd-tauStart)/4;

%The first one is an arbitrary integration constant and is just set to
%zero.
curN=0;
b(curN+1)=0;

%The first coeficient has a special formula.
if(n>=2)
    curN=1;
    b(curN+1)=normConst*(2*a(curN-1+1)-a(curN+1+1));

    %The recursion for most of the entries.
    for curN=2:(n-1)
        b(curN+1)=normConst*(a(curN-1+1)-a(curN+1+1))/curN;
    end
    %The formula for the second to last entry
    curN=n;
    b(curN+1)=normConst*(a(curN-1+1))/curN;
elseif(n==1)%If n=1.
    curN=n;
    b(curN+1)=normConst*(2*a(curN-1+1));
elseif(n==0)
    curN=n+1;
    b(curN+1)=normConst*(2*a(curN-1+1));
    return
end

%The coefficient for the last item has a special formula.
curN=n+1;
b(curN+1)=normConst*a(curN-1+1)/curN;

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
