function val=doubleFactorial(x)
%%DOUBLEFACTORIAL Evaluate x!!, which is the double factorial (also known
%                 as semifactorial) of x. If x is odd, then the double
%                 factorial is defined as x*(x-2)*(x-4)*...*5*3*1. If x is
%                 even, then it is x(x-2)*(x-4)*...*6*4*2. The double
%                 factorials of 0 and are defined to be 1.
%
%INPUTS:  x   A scalar or matrix, where the double factorial of each
%             element should be evaluated. The values must be positive
%             integers, or 0 or -1.   
%
%OUTPUTS: val The double factorial of each value of x.
%
%The gamma function/ factorial representations of the double factorial for
%even and odd values are taken from [1].
%
%The functions are implemented using logarithmic functions so as to
%minimize the effects of overflow with intermediate results.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Double Factorial." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/DoubleFactorial.html
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numEl=length(x(:));

val=zeros(size(x));
for curEl=1:numEl
    xCur=x(curEl);
    
    if(xCur==0||xCur==-1)
        val(curEl)=1;
    elseif(mod(xCur,2)==0)%If it is divisible by 2.
        k=xCur/2;
        val(curEl)=exp(k*log(2)+gammaln(k+1));
    else%Assume that it is odd.
        n=(xCur+1)/2;
        val(curEl)=exp(gammaln(xCur+2)-(n*log(2)+gammaln(n+1)));
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
