function aList=genAllMixRadGrayCode(maxVals)
%%GENALLMIXRADGRAYCODE This generates all tuples in a non-binary gray code.
%               Only one component is changed by +/-1 between successive
%               tuples in the output.
%
%INPUTS: maxVals An nX1 or 1Xn list of the maximum values of each element
%                of each tuple. Values in the ith element of the tuples go
%                from 0 to maxVals(i). n can be any length. An empty matrix
%                being passed leads to an empty matrix being returned.
%
%OUTPUTS: aList An nXnumTuples list of all of the tuples given in a gray
%               code ordering.
%
%This function implements Algorithm H from Section 7.2.1.1 of [1]. There
%are a total of numTuples=prod(maxVals+1) tuples.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(maxVals);
m=maxVals+1;

numTuples=prod(m);

aList=zeros(n,numTuples);

%Step H1, Initialize
a=zeros(n,1);
o=ones(n,1);
f=0:n;

for curTuple=1:numTuples
    %Step H2, visit the tuple.
    aList(:,curTuple)=a;

    %Step H3, Choose j.
    j=f(0+1);
    f(0+1)=0;

    %Step H4, Change j.
    if(j==n)
        break;
    end
    
    a(j+1)=a(j+1)+o(j+1);

    %Step H5
    if(a(j+1)==0||a(j+1)==m(j+1)-1)
        o(j+1)=-o(j+1);
        f(j+1)=f(j+1+1);
        f(j+1+1)=j+1;
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
