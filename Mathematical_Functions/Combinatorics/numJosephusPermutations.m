function val=numJosephusPermutations(n)
%%NUMJOSEPHUSPERMUTATIONS Obtain the number of unique length n Josephus
%    permutations. For n>=4, this is less then factorial(n).
%
%INPUTS: n The integer length of the Josephus permutation. n>=1.
%
%OUTPUTS: val A uint64 value holding the number of unique length n Josepus
%             permutations. For n>=42, this is clipped to 2^(64)-1,
%             because it is more than can be represented in a 64-bit
%             integer.
%
%The number of unique Josephus permutations is given in [1]. It is simply
%the least common multiple of the numbers 1:n.
%
%EXAMPLE:
%As noted in [1], for n=12, the number of Josephus permutations is 27,720.
%That is obtained here:
% val=numJosephusPermutations(12)
%
%REFERENCES
%[1] P. Schumer, "The Josephus problem: Once more around," Mathematics
%    Magazine, vol. 75, no. 1, pp. 12-17, Feb. 2002.
%
%February 2025 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n>=42)
    warning('Integer overflow occurred. The maximum 64 bit integer is returned.')
    val=uint64(18446744073709551615);%This is 2^64-1
    return;
end

val=lcm(uint64(1),uint64(2));
for k=3:n
    val=lcm(val,uint64(k));
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
