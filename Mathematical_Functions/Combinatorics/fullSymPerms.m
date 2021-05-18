function permVals=fullSymPerms(x)
%%FULLSYMPERMS This function generates all fully symmetric (multiset) 
%              permutations of the input sequence. That is, all possible
%              multiset permutations of the input sequence, with all
%              possible variations of sign (+/-) of the elements, making
%              sure not to repeat sequances, as might occur when dealing
%              with flipping the sign on 0. Such fully symmetric
%              permutation sets arise often when generating cubature
%              integration points, as in Chapter 8 of [1].
%
%INPUTS: x An NX1 vector.
%
%OUTPUTS: permVals An NXnumVals matrix of all unique fully symmetric
%                  permutations of the elements in x.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numNonzero=sum(x~=0);

%The special case of an all zero vector
if(numNonzero==0)
   permVals=x;
   return;
end

%Generate all multiset permutations.
multiSetPerms=genAllMultisetPermutations(x);

numMultiset=size(multiSetPerms,2);
xDim=size(x,1);

numCombos=2^numNonzero;

%Allocate space
permVals=zeros(xDim,numMultiset*numCombos);

curStart=1;
for curPerm=1:numMultiset
    permVals(:,curStart:(curStart+numCombos-1))=PMCombos(multiSetPerms(:,curPerm));
    curStart=curStart+numCombos;
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
