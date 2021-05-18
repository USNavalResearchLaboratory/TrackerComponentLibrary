function zList=genAllNestedZParenth(n,algorithm)
%%GENALLNESTEDZPARENTH0 Generate all nested sets of n parenthesis pairs,
%           recording only the indices of the left parantheses.
%
%INPUTS: n The number of nested paranthesis pairs; n>=0. Passing zero will
%          result in an empty matrix being returned.
% algorithm An optional parameter specifying the algorithm that should be
%          used to generate the parenthesis pairs. Possible values are:
%          0 (The default if omitted or an empty matrix is passed) Use
%            the algorithm of Problem 2 of Section 7.2.1.6 of [1].
%          1 Use algorithm N of Section 7.2.1.6 of [1].
%
%OUTPUTS: zList The nXnumPairs list of indices where the left parenthesis
%               goes.
%
%There is a total of CatalanNumber(n) possible arrangements of nested
%parentheses.
%
%EXAMPLE:
%To recreate the values in Table I of Section 7.2.1.6 of [1],
%one can run
% zList=genAllNestedZParenth(4)
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

if(n==0)
    zList=[];
    return;
end

switch(algorithm)
    case 0
        zList=genAllNestedZParenth0(n);
    case 1
        zList=genAllNestedZParenth1(n);
    otherwise
        error('Unknown algorithm specified.')
end

end

function zList=genAllNestedZParenth0(n)
%%GENALLNESTEDZPARENTH0 Generate all nested sets of n parenthesis pairs,
%           recording only the indices of the left parantheses using the
%           algorithm in Problem 2 of Section 7.2.1.6 of [1].
%
%REFERENCES
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

if(n==1)
    zList=1;
    return;
end

numBinTrees=CatalanNumber(n);
%Allocate space.
zList=zeros(n,numBinTrees);

%Step T1 Initialize
z=2*(0:n).'-1;

for curTree=1:numBinTrees
    %Step T2
    zList(:,curTree)=z(2:(n+1));
    
    %Step T3
    if(z(n-1+1)<z(n+1)-1)
        z(n+1)=z(n+1)-1; 
        continue;
    end
    
    %Step T4, find j.
    j=n-1;
    z(n+1)=2*n-1;
    while(z(j-1+1)==z(j+1)-1)
        z(j+1)=2*j-1;
        j=j-1;
    end
    
    if(j==1)
        break;
    end
    z(j+1)=z(j+1)-1;
end

end

function zList=genAllNestedZParenth1(n)
%%GENALLNESTEDZPARENTH0 Generate all nested sets of n parenthesis pairs,
%           recording only the indices of the left parantheses using
%           Algorithm N of Section 7.2.1.6 of [1].
%
%REFERENCES
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.

numBinTrees=CatalanNumber(n);
%Allocate space.
zList=zeros(n,numBinTrees);

%Step N1
z=2*(1:n).'-1;
g=2*(1:n).'-2;

for curTree=1:numBinTrees
    %Step N2
    zList(:,curTree)=z;
    j=n;

    %Step N3, Find j.
    while(z(j)==g(j))
        g(j)=bitxor(g(j),1);
        j=j-1;
    end

    %Step N4
    if(mod(g(j)-z(j),2)==0)
        z(j)=z(j)+2;
        continue;
    end
    
    %Step N5
    t=z(j)-2;
    if(t<0)
        return
    end

    if(t<=z(j-1))
        t=t+2*(t<z(j-1))+1;
    end
    z(j)=t;
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
