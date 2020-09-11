function y=polyValGen(p,x,direction)
%%POLYVALGEN Return the value of a scalar polynomial evaluated at the point
%       or points given in x. The polynomial can be specified in the format
%       y=p(1)*x^N+p(2)*x^(N-1)+...+p(N)*x+p(N+1)
%       or in the format
%       y=p(N+1)*x^N+p(N)*x^(N-1)+...+p(2)*x+p(1)
%       as specified by the direction option. This function is similar to
%       polyval, but allows one to choose between the formats without
%       calling flip to reverse the order of the elements in p.
%
%INPUTS: p An NX1 or 1XN set of real or complex polynomial coefficients;
%          wN>=1.
%        x A real or complex scalar point or a matrix of points at which
%          the scalar polynomial should be evaluated.
% direction An optional parameter specifying which of the two polynomial
%          formats is used. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) p(1) is
%            the coefficient of the highest order x term.
%          1 p(1) is the coefficient of the lowest order x term.
%
%OUTPUTS: y The value or values of the polynomial in p evaluated at the
%           point(s) in x. This is the same size as x, unless x is an empty
%           matrix, in which case y is 0.
%
%The function is evaluated using Horner's rule, which is described in [1].
%
%REFERENCES:
%[1] Weisstein, Eric W. "Horner's Rule." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/HornersRule.html
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(direction))
   direction=0; 
end

numP=length(p);

y=zeros(size(x));
if(direction==0)
    if(numP>0)
        y(:)=p(1);
    end
    
    for k=2:numP
        y=x.*y+p(k);
    end
else
    if(numP>0)
        y(:)=p(numP);
    end
    
    for k=(numP-1):-1:1
        y=x.*y+p(k);
    end
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
