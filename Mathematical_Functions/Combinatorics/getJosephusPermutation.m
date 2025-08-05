function vals=getJosephusPermutation(n,k)
%%GETJOSEPHUSPERMUTATION Obtain the kth Josephus permutation of length n.
% Given n items in a circle, one removed every kth item without
% replacement. The Josephus permutation is the sequence of the k items
% removed.
%
%INPUTS: n The total number of items.
%        k The kth item is taken each time. This is that number. Note that
%          k can be bigger than n, because it loops back.
%
%OUTPUTS: vals An nX1 vector holding values from 1 to n corrresponding to
%              the kth Josephus permutation of length n.
%
%The Josephus permutations are described in [1] and [2]. The algorithm is
%implied by solving Problem 14-2-b on page 318 of [3] making use of the
%concepts in the chapter. That is, it makes use of a priority search tree
%to speed things up. 
%
%EXAMPLE:
%This returns the permutation shown on page 16 of [1]. The  order is
%[15, 13, 11, 9, 7, 5, 3, 1, 14, 10, 6, 2, 12, 4, 8, 16]
% vals=getJosephusPermutation(16,720719)
%
%REFERENCES
%[1] P. Schumer, "The Josephus problem: Once more around," Mathematics
%    Magazine, vol. 75, no. 1, pp. 12-17, Feb. 2002.
%[2] Weisstein, Eric W. "Josephus Problem." From MathWorld--A Wolfram Web
%    Resource. https://mathworld.wolfram.com/JosephusProblem.html
%[3] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein,
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%
%February 2025 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

theTree=AVLTree();
for i=1:n
    theTree.insert(KeyVal(i,[]));
end

vals=zeros(n,1);
pos=0;
numAdded=0;
k=k-1;
while(theTree.numInTree>1)
    pos=mod(pos+k,theTree.numInTree);
    numAdded=numAdded+1;
    vals(numAdded)=theTree.findIthSmallestElement(pos+1).key;
    theTree.remove(vals(numAdded));
end
vals(n)=theTree.keyVal.key;

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
