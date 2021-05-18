function val=EulerianNumber(n,k)
%%EULERIANNUMBER Compute an Eulerian number. This is the number of
%                permutation of 1,2,...n that have k permutation ascents.
%                Given a permutation a_1...a_n, a permutation asent is if
%                a_i<a_i+1.
%
%INPUTS: n The integer length of the permutation. n>=0.
%        k The integer number of ascents desired. k>=0.
%
%OUTPUTS: val The number of permutations of length n having k ascents.
%
%The Eulerian numbers are computed using the recursion given in [1]. THe
%algorithm is slow for large values of n. In the manner that Binomial
%coefficients form Pascal's triangle with each row summing to 2^(row
%number), Eulerian coefficients form Euler's triangle with each row summing
%to (row number)!.
%
%REFERENCES:
%[1] Weisstein, Eric W. "Eulerian Number." From MathWorld--A Wolfram Web
%    Resource. http://mathworld.wolfram.com/EulerianNumber.html
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
if(k<=0||k>n)
    val=0;
    return;
end
    
prevRow=zeros(k+1,1);
curRow=zeros(k+1,1);

prevRow(1)=1;

for nCur=2:n
    curRow(1)=1;
    for kCur=2:(nCur-1)
        curRow(kCur)=(nCur-kCur+1)*prevRow(kCur-1)+kCur*prevRow(kCur);
    end
    curRow(nCur)=1;
    prevRow=curRow;
end

val=curRow(k);
    
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
