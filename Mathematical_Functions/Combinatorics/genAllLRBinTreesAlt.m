function [lList,rList]=genAllLRBinTreesAlt(n)
%%GENALLLRBINTREESALT Generate representations of all binary trees having n
%                     nodes. l(0+1) gives the root node index. l(i+1) and
%                     r(i) then give the left and right descendents of node
%                     i. If there is a zero, then there is no descendent.
%                     The structure of the (left-right) relationships is
%                     what is unique, not the specific nodes. That is, the
%                     identity of the root node does not really matter. For
%                     example, node 1 will never have a left neighbor and
%                     node 4 will never have a right neighbor. However,
%                     independent of node labels, all binary trees are
%                     formed. The format of the output differs from the
%                     function genAllLRBinTrees.
%
%INPUTS: n The number of nodes that will be in the trees. n>=1.
%
%OUTPUTS: lList The (n+1)XnumTrees lists of left links, where lList(i+1,k)
%               is the left link of node i in tree k with lList(0+1,k)
%               being the root.
%         rList The nXnumTrees lists of right links.
%
%This function implements Algorithm L of Section 7.2.1.6 of [1]. The number
%of binary trees is CatalanNumber(n).
%
%EXAMPLE:
%To recreate the list of links in Table 3 of Section 7.2.1.6 of [1], one
%can use
% [lList,rList]=genAllLRBinTreesAlt(4)
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numBinTrees=CatalanNumber(n);

%Allocate space
lList=zeros(n+1,numBinTrees);
rList=zeros(n,numBinTrees);
kList=zeros(n,numBinTrees);

k=zeros(n,1);
o=zeros(n+1,1);
l=zeros(n+1,1);
r=zeros(n,1);

%Step L1
l((1:n)+1)=0;
r(1:(n-1))=2:n;
k(1:n)=0:(n-1);
o((1:n)+1)=-1;
l(0+1)=1;
o(0+1)=1;
r(n)=0;

for curTree=1:numBinTrees
    %Step L2
    lList(:,curTree)=l;
    rList(:,curTree)=r;
    kList(:,curTree)=k;
    j=n;
    p=0;    
    
    while(1)
        %Step L3
        if(o(j+1)>0)
            m=l(j+1);

            if(m~=0)
               %Step L5 
                if(j==0)
                    return;
                end
                l(j+1)=r(m);
                r(m)=j;
                k(j)=m;
                x=k(m);
                if(x==0)
                    l(p+1)=m;
                else
                    r(x)=m;
                end
                break;
            else
                %p=j;
                o(j+1)=-o(j+1);
                j=j-1;
                continue;
            end
        else%o(j+1)<0
            m=k(j);

            if(m~=0)
                %Step L4
                r(m)=l(j+1);
                l(j+1)=m;
                x=k(m);
                k(j)=x;
                if(x==0)
                   l(p+1)=j;
                else
                   r(x)=j;
                end
                break;
            else
                p=j;
                o(j+1)=-o(j+1);
                j=j-1;
                continue;
            end
        end
    end
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
