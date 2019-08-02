function theRank=rankCombination(combo,firstElMostSig,startVal)
%%RANKCOMBINATION Obtain the colexicographic order (the rank) of a
%                 particular combination starting from 0, where the first
%                 element of combo is the least significant, unless
%                 otherwise specified.
%
%INPUTS: combo A vector where each elements holds a number of the
%              combination in INCREASING order unless firstElMostSig is
%              true. The numbers are positive integers starting from
%              startVal. For example, [2;1;0] and [5;4;3] are valid
%              combinations, but [3;4;5] is an invalid combination (unless
%              firstElMostSig is true), because the elements are not in
%              decreasing order.
% firstElMostSig An optional parameter specifying whether the first
%              element is the most or least significant. The default if
%              omitted or an empty matrix is passed is false (the first
%              element is the least significant).
%     startVal This is zero or 1, indicating which value the value at which
%              the elements in combo can start. The default if omitted or
%              an empty matrix is passed is 0.
%
%OUTPUTS: theRank The rank of the combination in a colexicographic ordering
%                 of combinations, counting from 0.
%
%The algorithm is an implementation of the combinatorial number system of
%degree m, where m is the length of the vector combo. The system is
%described in Chapter 7.2.1.3 of [1].
%
%REFERENCES:
%[1] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(isempty(combo))
   theRank=[];
   return;
end

if(nargin<3||isempty(startVal))
    startVal=0;
end

if(nargin<2||isempty(firstElMostSig))
    firstElMostSig=false;
end

combo=combo-startVal;

m=length(combo);
if(firstElMostSig==false)
    combo=flipud(combo(:));
end

theRank=0;
for curIdx=1:m
    theRank=theRank+binomial(combo(curIdx),m-curIdx+1);
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
