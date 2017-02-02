function val=FibonacciNum(n)
%%FIBONACCINUM  Determine the nth Fibonacci number. Fibonacci numbers are
%               the sequence 1, 1, 2, 3, 5, 8, ... for n=1,2,3,4,5,6,...
%               where F_{n}=F_{n-1}+F_{n-2} determines the next value. A
%               non-recursion expression is used here so that large
%               Fibonacci numbers can be quickly found.
%
%INPUTS: n   A scalar or matrix of integer Fibonacci number positions, n>0.
%
%OUTPUTS: val The values of the Fibonacci numbers at the positions given in
%             n. These should be exact for n<70.
%
%The recursive and non-recursive formulae are taken from [1]. The recursive
%formula is used up to n=78, as after that, finite precision errors
%dominate. For n=79 onewards, the non-recusive formula is used. The formula
%is written in a slightly more convoluted manner than one might expect
%using log and exp to try to minimize finite precision errors as
%much as possible. The issue comes with the exponent of phi. The
%non-recusive formula is not used when the recursive formula is exact as
%the non-recusive formula loses a few bits of precision before the recusive
%formula (in the n=70's.) However, the non-recusive formula is faster and
%has more bits of precision for very large values of n.
%
%REFERENCES:
%[1] Chandra, Pravin and Weisstein, Eric W. "Fibonacci Number." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/FibonacciNumber.html
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Allocate space
val=zeros(size(n));

if(any(n~=fix(n)|imag(n)~=0|n<=0))
   error('n has the wrong format')
end

%Set the starting cases.
val(n==1|n==2)=1;

%Select values that are small enough so that the recursion can be used
%without a loss of precision
sel=(n<=78)&(n>2);

y=n(sel);
%Use the definition involving the recursion.
if(~isempty(y))
    numY=numel(y);
    res=zeros(numY,1);
    
    %We sort it in increasing order, so that we do not have to go back.
    y=sort(y);
    
    fCur=1;
    fPrev=1;
    curY=1;
    curN=3;
    while(curY<=numY)
        f=fCur+fPrev;
        
        if(y(curY)==curN)
            res(curY)=f;
            curY=curY+1;
        end
        
        fPrev=fCur;
        fCur=f;
        
        curN=curN+1;
    end
    
    val(sel)=res;
end
y=n(~sel);
%Deal with big numbers
if(~isempty(y))
    %The golden ratio.
    phi=(1+sqrt(5))/2;

    val(~sel)=round(exp(y*log(phi)-log(sqrt(5))));
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
