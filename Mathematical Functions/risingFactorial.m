function val=risingFactorial(x,n,algorithm)
%%RISINGFACTORIAL Evaluate the rising factorial <x>_n, also denoted x^(n),
%                 and also known as the Pochhammer symbol. The rising
%                 factorial <x>_n=x*(x+1)*(x+2)*...*(x+n-1), where <x>_n=1.
%                 It can be expressed as gamma(x+n)/gamma(x), which makes
%                 it continuous. Other definitions can arise (but are not
%                 offered here) when considering negative non-integer
%                 values of n.
%
%OUTPUTS: x A scalar or matrix of first values in the product series of
%           x*(x+1)*...*(x+n-1). x must be real.
%         n The scalar or matrix n value in the product series of
%           x*(x+1)*...*(x+n-1). n must be real.
% algorithm Optionally, one can specify the algorithm to use. If this
%           parameter is omitted, then algorithm 0 is used for positive
%           integer values of n, then algorithm 1 if x>0 and x+n>0,
%           otherwise algorithms 2 is used. The possible values are:
%           0 Use the product identity that is given in [1].
%           1 Use the expression in terms of gamma functions from [1], but
%             use gammaln to reduce the effects of possible overflows.
%           2 Use the expression in terms of gamma functions from [1]
%             directly using the gamma function.
%
%OUTPUTS: val The rising factorial <x>_n for each x and the corresponding
%             n or for each n and a single x depending on what is a matrix
%             and what is scalar.
%
%Rising factorials often arise when dealing with combinatorial problems and
%hypergeometric functions.
%
%The rising factorial is discussed in [1].
%Note that for positive integer n, risingFactorial(1,n)=factorial(n)
%
%REFERENCES:
%[1] Weisstein, Eric W. "Rising Factorial." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/RisingFactorial.html
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=zeros(size(x));
numVals=length(x(:));

for curVal=1:numVals
    xCur=x(curVal);
    
    if(nargin<3||isempty(algorithm))
        if(fix(n)==n&&n>=0)
            algorithm=0;
        elseif(xCur>0&&xCur+n>0)
            algorithm=1;
        else
            algorithm=2;
        end
    end

    switch(algorithm)
        case 0%The product identity for positive integer n
            val(curVal)=1;
            for k=1:n
                val(curVal)=(xCur+k-1)*val(curVal);
            end
        case 1
            val(curVal)=exp(gammaln(xCur+n)-gammaln(xCur));
        case 2
            val(curVal)=gamma(xCur+n)/gamma(xCur);
        otherwise
            error('Unknown algorithm specified');
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
