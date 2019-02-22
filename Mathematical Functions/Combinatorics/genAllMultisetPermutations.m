function thePerms=genAllMultisetPermutations(a,algorithm)
%%GENALLMULTISETPERMUATIONS Generate all unique permuations of the elements
%               in a multiset. A multiset is a set of elements where some
%               elements can be repeated. The algorithm chosen affects the
%               ordering of the outputs.
%
%INPUTS:  a An nX1 or 1Xn vector of the elements in the multiset, including
%           repetitions. The elements of a must be finite; a cannot contain
%           any NaNs.
% algorithm An optional parameter specifying the algorithm to use. possible
%           values are:
%           0 (The default if this parameter is omitted or an empty matrix
%             is passed). Use algorithm L of Chapter 7.2.1 of [1]. This
%             produces the multiset permutations in lexicographic order.
%             This algorithm comes from Narayana Pandita in 14th century
%             India.
%           1 Use the algorithm of [2], implemented with arrays rather than
%             explicit data structures for linked lists. The permutations
%             of the multiset are ordered such that neighboring
%             permutations are related by a minimal change. This is not
%             lexicographic order.
%
%OUTPUS: thePerms A nXnumPerm matrix of the permutations of the multiset.
%                 The ordering is specified by the algorithm used.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%[2] A. Williams, "Loopless generation of multiset permutations by prefix
%    shifts," in Proceedings of the Twentieth Annual ACM-SIAM Symposium
%    on Discrete Algorithms, New York, NY, 4-6 Jan. 2009, pp. 987-996.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
   algorithm=0; 
end

switch(algorithm)
    case 0
        thePerms=genAllMultisetPermsAlgL(a);
    case 1
        thePerms=multisetPermsMinOrder(a);
    otherwise
        error('Unknown algorithm specified');
end
end

function thePerms=multisetPermsMinOrder(E)
%%GENALLMULTISETPERMUTATIONS Generate all unique permutations of the
%               elements in a multiset using the minimum order algorithm of
%               [1].
%
%INPUTS: E An nX1 or 1Xn vector of the elements in the multiset, including
%          repetitions. The elements of E must be finite; E cannot contain
%          any NaNs.
%
%OUTPUS: thePerms A nXnumPerm matrix of the permutations of the multiset
%              ordered such that neighboring permutations are related by a
%              minimal change. This is not lexicographic order. 
%
%This function implements Algorithm 1 in [1]. However, this implementation
%uses arrays rather than explicit data structures for the linked list. The
%function numMultisetPermutations provides the number of multiset
%permutations for a given vector of elements in the multiset.
%
%REFERENCES:
%[1] A. Williams, "Loopless generation of multiset permutations by prefix
%    shifts," in Proceedings of the Twentieth Annual ACM-SIAM Symposium
%    on Discrete Algorithms, New York, NY, 4-6 Jan. 2009, pp. 987-996.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(E);

EList=sort(E,'descend');
EList=EList(:);

numPerms=numMultisetPermutations(EList);
thePerms=zeros(n,numPerms);

nextList=[2:n,0];
head=1;
i=n-1;
afterI=n;

%Visit the first value, which is the list sorted in descending order.
thePerms(:,1)=EList;

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
        thePerms(curIdx,curPerm)=EList(curNode);
        curNode=nextList(curNode);
    end
    curPerm=curPerm+1;
end

end

function thePerms=genAllMultisetPermsAlgL(a)
%%GENALLMULTISETPERMSALGL Generate all multiset permutations in
%               lexicographic order using algorithm L of Chapter 7.2.1.2 of
%               [1].
%
%INPUTS: a An nX1 or 1Xn vector of the elements in the multiset, including
%          repetitions. The elements of a must be finite; a cannot contain
%          any NaNs.
%
%OUTPUS: thePerms A nXnumPerm matrix of the permutations of the multiset
%                 given in lexicographic order.
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=length(a);

numPerms=numMultisetPermutations(a);
thePerms=zeros(n,numPerms);

%The first element is the a(0) term added for convenience.
a=sort(a(:),'ascend');

curPerm=1;
while(1)
    %Step L1, visit the permutation
    thePerms(:,curPerm)=a(1:n);
    if(curPerm==numPerms)
        break;
    end

    %Step L2, find j. Here, we do not need a check for j=0, because we
    %previously checked for termination.
    j=n-1;
    while(a(j)>=a(j+1))
        j=j-1;
    end

    %Step L3
    l=n;
    while(a(j)>=a(l))
        l=l-1;
    end
    temp=a(j);
    a(j)=a(l);
    a(l)=temp;
    
    %Step L4
    a((j+1):n)=a(n:-1:(j+1));

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
