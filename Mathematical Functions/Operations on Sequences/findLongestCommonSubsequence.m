function xCommon=findLongestCommonSubsequence(x,y)
%%FINDLONGESTCOMMONSUBSEQUENCE Given two vectors of numbers, characters,
%              etc., find the longest sequence of elements in both of them.
%              The elements must appear in order need NOT BE CONSECUTIVE.
%
%INPUTS:      x,y  Two vectors.
%
%OUTPUTS: xCommon  The longest common subsequence shared by the vectors.
%                  The elements in the subsequence need not be consecutive.
%                  For example, the longest common subsequence of "dunkin"
%                  and "doughnuts" is "dun".
%
%The algorithm is taken from Chapter 15.4 of [1], but modified so that the
%sequence comes out in order and to deal with Matlab addressing things from
%1 instead of 0.
%
%REFERENCES:
%[1] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein,
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=length(x);
n=length(y);

c=zeros(m+1,n+1);
for i=m:-1:1
    for j=n:-1:1
        if(x(i)==y(j))
            c(i,j)=c(i+1,j+1)+1;
        else
            c(i,j)=max(c(i,j+1),c(i+1,j));
        end
    end
end

i=1;
j=1;
%Setting xCommon to x rather than alocating with zeros assured that xCommon
%will have the same type and shape as x. This is good when x is a character
%string.
xCommon=x;
numFound=0;
while(i<=m&&j<=n)
    if(x(i)==y(j))
       numFound=numFound+1;
       xCommon(numFound)=x(i);
       i=i+1;
       j=j+1;
    elseif(c(i+1,j)>=c(i,j+1))
        i=i+1;
    else
        j=j+1; 
    end
end

%Get rid of extra elements in xCommon.
xCommon=xCommon(1:numFound);
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
