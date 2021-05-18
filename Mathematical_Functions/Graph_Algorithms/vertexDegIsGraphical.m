function boolVal=vertexDegIsGraphical(d)
%%VERTEXDEGISGRAPHICAL Given a list of the degrees of vertices in a desired
%           non-dicrectional graph, determine whether a real graph exists
%           (a set of edges exists) that can satisfy the degree profile.
%            This does not explicitly find the graph.
%
%INPUTS: d An nX1 or 1Xn vector of the degrees of the n vertices. These are
%          integer values >=0. They represent the number of edges touching
%          the nodes.
%
%OUTPUTS: boolVal This is true if the degree list is graphical and false
%                 otherwise.
%
%This function implements Problem 105 of Chapter 7 of [1].
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(d);
d=sort([d(:);0],'descend');

if(mod(sum(d),2)==1)
    boolVal=false;
    return;
end

c=zeros(d(1),1);

%Step H1
k=d(1);
j=0;
while(k>0)
    j=j+1;
    while(k>d(j+1))
        c(k)=j;
        k=k-1;
    end
end
%If all of the d's are zero.
if(j==0)
    boolVal=true;
    return;
end

s=0;
while(d(s+1)>=s+1&&s+1<=min(n,d(1)))
    s=s+1;
end

if(s==0)
    boolVal=false;
    return;
end

dSum=0;
cSum=0;
for k=1:s
    dSum=dSum+d(k);
    cSum=cSum+c(k);
    
    if(dSum>cSum-k)
        boolVal=false;
        return;
    end
end

boolVal=true;

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
