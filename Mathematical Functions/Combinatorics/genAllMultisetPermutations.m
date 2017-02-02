function perms=genAllMultisetPermutations(E)
%%GENALLMULTISETPERMUTATIONS Generate all unique permutations of the
%               elements in a multiset. A multiset is a set of elements
%               where some elements are repeated.
%
%INPUTS: E An nX1 or 1Xn vector of the elements in the multiset, including
%          repetitions. The elements of E must be finite; E cannot contain
%          any NaNs.
%
%OUTPUS: perms A nXnumPerm matrix of the permutations of the multiset
%              ordered such that neighboring permutations are related by a
%              minimal change. This is not lexicographic order. 
%
%This function implements Algorithm 1 in [1]. However, this implementation
%uses arrays rather than explicit data structures for the linked list. The
%function numMultisetPermutations prodives the number of multiset
%permutations for a fiven vector of elements in the multiset.
%
%REFERENCES:
%[1] A. Williams, "Loopless generation of multiset permutations by prefix
%    shifts," in Proceedings of the Twentieth Annual ACM-SIAM Symposium
%    on Discrete Algorithms, New York, NY, 4?6 Jan. 2009, pp. 987-996.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(E);

EList=sort(E,'descend');
EList=EList(:);

numPerms=numMultisetPermutations(EList);
perms=zeros(n,numPerms);

nextList=[2:n,0];
head=1;
i=n-1;
afterI=n;

%Visit the first value, which is the list sorted in descending order.
perms(:,1)=EList;

curPerm=2;
while(nextList(afterI)~=0||EList(afterI)<EList(head))
    if(nextList(afterI)~=0&&EList(i)>=EList(nextList(afterI)))
        beforeK=afterI;
    else
        beforeK=i;
    end
    
    k=nextList(beforeK);
    nextList(beforeK)=nextList(k);
    nextList(k)=head;
    if(EList(k)<EList(head))
       i=k; 
    end
    afterI=nextList(i);
    head=k;

    %Visit the permutation. This means starting at head and going until
    %reaching the end.
    curNode=head;
    for curIdx=1:n
        perms(curIdx,curPerm)=EList(curNode);
        curNode=nextList(curNode);
    end
    curPerm=curPerm+1;
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
