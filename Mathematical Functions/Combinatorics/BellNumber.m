function val=BellNumber(n)
%%BELLNUMBER Compute the nth Bell number. This is the number of ways a set
%            of n elements can be partitioned into nonempty subsets. The
%            number is computed using a recursion that will not overflow
%            unless the actual number overflows. When n is very large, the
%            function can be slow.
%
%INPUTS:  n The integer total number of items in the set n>=0.
%
%OUTPUTS: val The number of ways of partitioning a set of n items into
%             non-empty subsets.
%
%Bell numbers are defined in terms of Stirling numbers of the second king
%in [1]. This implementation uses the recurrent relation for Stirling
%numbers of the second kind from [2], summing over the final recurrence
%row.
%
%The more commonly seen expression in terms of a sum of binomials is not
%desirable, because overflows will occur even when the final number does
%not overflow.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Bell Number." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/BellNumber.html
%[2] Weisstein, Eric W. "Stirling Number of the Second Kind." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/StirlingNumberoftheSecondKind.html
%
%January 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

curN=zeros(n+1,1);
prevN=zeros(n+1,1);

prevN(0+1)=1;%{0,0}=1
for nCur=1:n
    for kCur=1:nCur
        curN(kCur+1)=prevN(kCur-1+1)+kCur*prevN(kCur+1);
    end
    prevN=curN;
end
val=sum(prevN);

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
