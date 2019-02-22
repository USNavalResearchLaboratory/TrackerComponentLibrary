function val=fallingFactorial(x,n)
%%FALLINGFACTORIAL Evaluate the falling factorial (x)_n, which for n>=0 is
%                  also known as the binomial polynomial, the lower
%                  factorial, the factorial power, and the falling
%                  factorial power. It is defined
%                  (x)_n=x*(x-1)*(x-2)*...*(x-(n-1)). It can be expressed
%                  as gamma(x+1)/gamma(x-n+1), which makes it continuous.
%                  It is implemented using logarthmic functions to reduce
%                  problems with overflows in intermediate results.
%
%INPUTS: x A scalar or matrix of first values in the product series of
%          x*(x-1)*...*(x-(n-1)). x cannot be negative.
%        n The scalar or matrix n value in the product series of
%          x*(x-1)*...*(x-(n-1)). If x is a matrix, n can be a scalar or a
%          matrix of the same size. If x is a scalar, n can be a scalar or
%          a matrix.
%
%OUTPUTS: val The falling factorial (x)_n for each x and the corresponding
%             n/ for each n and a single x depending on what is a matrix
%             and what is scalar. Note that if n>x, the falling factorial
%             is defined to be 0.
%
%The falling factorial is discussed in [1]. Note that
%fallingFactorial(n,n)=factorial(n);
%
%REFERENCES:
%[1] Weisstein, Eric W. "Falling Factorial." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/FallingFactorial.html
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=exp(gammaln(x+1)-gammaln(max(0,x-n+1)));

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
