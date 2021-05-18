function [vals,coeffs]=acosDerivVal(z,n,coeffs)
%%ACOSDERIVVAL Return the value of the nth derivative of acos(x) with
%             respect to x  evaluated at x=z for real or imaginary z.
%
%INPUTS: z A matrix of real or complex values at which the nth derivative
%          of the arccosine function is desired. If an empty matrix is
%          passed, then just coeffs will be returned.
%        n The number of derivatives to take. n>=0.
%   coeffs The computation of the derivatives involves determining the
%          structure of a polynomial. If this function has been run before
%          for a given n value, then the polynomial can be passed back and
%          is not determined again. Otherwise, this can be omitted or an
%          empty matrix can be passed. This only makes a difference for
%          n>=12.
%
%OUTPUTS: vals The value of the nth derivative of the arccosine function
%              taken at all of the points in z. vals has the same
%              dimensions as z.
%       coeffs A vector of polynomial coefficients that can be passed back
%              to this function when evaluating with the same n to speed it
%              up. 
%
%For n=0, vals=acos(z). Subsequent derivatives are computed using the
%relation that the derivative of the arccosine function is the negative of
%the derivative of the arcsine function, so we can just call asinDerivVal.
%
%May 2018 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n==0)
    vals=acos(z);
    coeffs=[];
    return;
end

if(nargin<3)
    coeffs=[];
end

[vals,coeffs]=asinDerivVal(z,n,coeffs);
vals=-vals;
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
