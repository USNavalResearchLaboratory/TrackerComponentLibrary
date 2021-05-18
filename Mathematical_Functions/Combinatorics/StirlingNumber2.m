function val=StirlingNumber2(n,k)
%%STIRLINGNUMBER2 Compute the Stirling number of the second kind {n,k},
%                 which is also called a Stirling set number. This is the
%                 number of ways to partition a set of n objects into k
%                 non-empty subsets. The number is computed using a
%                 recursion that will not overflow unless the actual number
%                 overflows. When n is very large, the function can be
%                 slow.
%                 
%INPUTS: n The integer total number of items in the set n>=0.
%        k The desired integer number of non-empty subsets of the set.
%          k>=0.
%
%OUTPUTS: val The number of ways of partitioning n items into k subsets. If
%             an overflow occurs, this will be infinite.
%
%The implementation uses the recurrence relation from [1].
%
%The more commonly seen expression in terms of a sum of binomials is not
%desirable, because overflows will occur even when the final number does
%not overflow.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Stirling Number of the Second Kind." From
%    MathWorld--A Wolfram Web Resource.
%    http://mathworld.wolfram.com/StirlingNumberoftheSecondKind.html
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n<k)
    val=0;
    return;
end

curN=zeros(k+1,1);
prevN=zeros(k+1,1);

prevN(0+1)=1;%{0,0}=1
for nCur=1:n
    for kCur=1:min(nCur,k)
        curN(kCur+1)=prevN(kCur-1+1)+kCur*prevN(kCur+1);
    end
    prevN=curN;
end
val=prevN(k+1);

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
