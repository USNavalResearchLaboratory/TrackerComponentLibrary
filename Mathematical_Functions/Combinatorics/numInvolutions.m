function val=numInvolutions(n)
%%NUMINVOLUTIONS Determine the number of involutions of length n. An
%                involution is a permution of items 1:n such that the
%                permutation is its own inverse. The number of involutions
%                also equals the number of possible Young tableau having n
%                cells.
%
%INPUTS: n A positive integer >=0.
%
%OUTPUTS: val The number of involutions of length n. For n=0, this is
%             defined as 1.
%
%This function uses the recurrence relation given in Equation 40 in Section
%5.1.4 of [1].
%
%REFERENCES:
%[1] D. Knuth, The Art of Computer Programming: Sorting and Searching, 2nd
%    ed. Boston: Addison-Wesley, 2018, vol. 3.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n<0)
    error('n must be >=0.')
end

if(n<=1)
    val=1;
    return;
end

an1=1;
an=1;
for k=2:n
    an2=an1;
    an1=an;
    an=an1+(k-1)*an2; 
end
val=an;
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
