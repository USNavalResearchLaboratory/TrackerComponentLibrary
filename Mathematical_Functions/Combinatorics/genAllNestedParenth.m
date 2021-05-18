function pStrings=genAllNestedParenth(n)
%%GENALLNESTEDPARENTH Generate all possible strings (character arrays) of
%              nested parentheses where there are n pairs of parentheses.
%              The possibilities are given in lexicographic order, where
%              ')' is considered lexicographically smaller than ')'. These
%              correspond to all possible forests of binary trees with n
%              nodes.
%
%INPUTS: n The number of nodes in the tree (The number of parenthesis
%          pairs); n>=1.
%
%OUTPUTS: pStrings A (2*n)XnumTrees character array of all possible valid
%                  parenthesis arrangements.
%
%This function implements Algorithm P in Section 7.2.1.6 of [1]. The case
%of n==1 is handles as a special case. A possible arrangement of
%parentheses is only valid if an opening ( always comes before its closing
%). Note that numTrees=CatalanNumber(n).
%
%EXAMPLE:
%Here, we generate all of the parenthesis pairs that are given in Table 1
%of Chapter 7.2.1.6 of [1].
% pStrings=genAllNestedParenth(4).'
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(n==1)
   pStrings="()".';
   return;
end

numBinTrees=CatalanNumber(n);

%Allocate space for the results.
pStrings=repmat(' ',2*n,numBinTrees);

%The added one if for an a(0). Here, we allocate space.
a=repmat(' ',2*n+1,1);

%Step P1, initialization.
%The a(0) entry, but indexation is from 1 here.
a(0+1)=')';
for k=1:n
   a(2*k-1+1)='(';
   a(2*k+1)=')';
end
m=2*n-1;

for curString=1:numBinTrees
    %Step P2, visit the string.
    pStrings(:,curString)=a(2:(2*n+1));

    %Step P3
    a(m+1)=')';
    if(a(m-1+1)==')')
        a(m-1+1)='(';
        m=m-1;
        continue;
    end
    
    %Step P4, find j.
    j=m-1;
    k=2*n-1;
    while(a(j+1)=='(')
        a(j+1)=')';
        a(k+1)='(';
        j=j-1;
        k=k-2;
    end

    %Step P5
    if(j==0)
        break
    end
    
    a(j+1)='(';
    m=2*n-1; 
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
