function [permList,invPermList]=genAllTopoConstPerms(n,orderReqMat)
%%GENALLTOPOCOSTPERMS Generate all permutations that obey certain
%           topological sorting restrictions. The restructions allowed are
%           to maintain the original ordering of some elements. If
%           orderReqMat(i,j)=1 (with i<j), then  it is required that the
%           value i appears before the value j in the permutations. Inverse
%           permutations are simultaneously generated and in such
%           instances, if a is an permutation vector, then the contrains
%           ensure that a(i)<a(j). Note that it is not possible to make an
%           order constraint orderReqMat(i,j)=1 for j<i. Thus, the initial
%           permutation of 1:n is always feasible.
%
%INPUTS: n The number of items in the permutations.
% ordReqMat An nXn matrix indicating the constraints. Only values of
%           orderReqMat(i,j) for i<j are considered. If orderReqMat(i,j)=1,
%           then the value i will always appear before j in the
%           permutations.
%
%OUTPUTS: permList The nXnumPerm set of order-constrained permutations.
%      invPermList The nXnumperm set of inverse permutations of those in
%                  permList.
%
%This function implements Algorithm V in Section 7.2.1.2 of [1].
%
%EXAMPLE:
%Here, we implement the example Given in Section 7.2.1.2 of [1], where all
%42 Young tableau of size 3X3 are generated from constrained inverse
%permutations.
% n=9;
% orderReqMat=zeros(9,9);
% orderReqMat(1,2)=1;
% orderReqMat(2,3)=1;
% orderReqMat(4,5)=1;
% orderReqMat(5,6)=1;
% orderReqMat(7,8)=1;
% orderReqMat(8,9)=1;
% orderReqMat(1,4)=1;
% orderReqMat(4,7)=1;
% orderReqMat(2,5)=1;
% orderReqMat(5,8)=1;
% orderReqMat(3,6)=1;
% orderReqMat(6,9)=1;
% [~,invPermList]=genAllTopoConstPerms(n,orderReqMat)
%Each of the permutations in invPermList corresponds to a Young tablea
%shown in [1].
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4A: Combinatorial
%    Algorithms, Part I, Boston: Addison-Wesley, 2011.
%
%October 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%At most, these will hold factorial(n) items. However, the total can be
%much smaller. We thus only allocate space for at most 5000 elements.
%Matlab can expand the arrays, if necessary.
maxAlloc=min(factorial(n),5000);
permList=zeros(n,maxAlloc);
invPermList=zeros(n,maxAlloc);

%Step V1
a=0:n;
aPrime=0:n;

curPerm=1;
while(1)
    %Step V2, visit the permutation.
    permList(:,curPerm)=a(2:(n+1));
    invPermList(:,curPerm)=aPrime(2:(n+1));
    curPerm=curPerm+1;
    k=n;

    while(1)
        %Step V3, can k move left?
        j=aPrime(k+1);
        l=a(j-1+1);

        if(~(l==0||orderReqMat(l,k)))
            %Step V4, Move it.
            a(j-1+1)=k;
            a(j+1)=l;
            aPrime(k+1)=j-1;
            aPrime(l+1)=j;
            break;%Go to Step V2
        end

        %Step V5, put k back.
        while(j<k)
            l=a(j+1+1);
            a(j+1)=l;
            aPrime(l+1)=j;
            j=j+1;
        end
        a(k+1)=k;
        aPrime(k+1)=k;
        k=k-1;
        if(k<=0)
            %Shrink to fit the actual number of permutations found.
            permList=permList(:,1:(curPerm-1));
            invPermList=invPermList(:,1:(curPerm-1));
            return;
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
