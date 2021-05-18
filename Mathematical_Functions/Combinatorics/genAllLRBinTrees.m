function [lList,rList]=genAllLRBinTrees(n)
%%GENALLLRBINTREES Generate representations of all binary trees having n
%                  nodes. The first node is indexed as 1 and is not
%                  explicitly given in the output. The output lists which
%                  node branches left or right from the previous node. The
%                  nodes are not ordered. Left and right node branches are
%                  given if l is the left branches and r is the right
%                  branches, then l(k) and r(k) are the nodes on the left
%                  and right of node k. If no node is on a branch, then
%                  l(k) and/or r(k) will be 0. The branch lists for all
%                  such trees are given in lList and rList.
%
%INPUTS: n The number of nodes that will be in the trees. n>=1.
%
%OUTPUTS: lList The nXnumTrees lists of left links. The final element in
%               each list is zero.
%         rList The nXnumTrees lists of right links. The final element in
%               each list is zero.
%
%This function implements Algorithm B of Section 7.2.1.6 of [1]. The number
%of binary trees is CatalanNumber(n).
%
%EXAMPLE:
%To recreate the list of links in Table 2 of Section 7.2.1.6 of [1], one
%can use
% [lList,rList]=genAllLRBinTrees(4)
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numBinTrees=CatalanNumber(n);

lList=zeros(n,numBinTrees);
rList=zeros(n,numBinTrees);

%Step B1, Initialize
l=[2:n,0,1].';%Length n+1
r=zeros(n,1);

for curTree=1:numBinTrees
    %Step B2
    lList(:,curTree)=l(1:n);
    rList(:,curTree)=r;
    
    %Step B3, find j.
    j=1;
    while(l(j)==0)
        r(j)=0;
        l(j)=j+1;
        j=j+1;
        if(j>n)
            return;
        end
    end

    %Step B4, find k and y.
    y=l(j);
    k=0;
    while(r(y)>0)
       k=y;
       y=r(y);
    end
    
    %Step B5, promote y.
    if(k>0)
        r(k)=0; 
    else
        l(j)=0;        
    end
    r(y)=r(j);
    r(j)=y;
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
