function theRank=rankCombination(combo,firstElMostSig,startVal,algorithm,n)
%%RANKCOMBINATION Obtain the position in a sequence (the rank) of a
%                 particular combination starting from 0, where the first
%                 element of combo is the least significant, unless
%                 otherwise specified. The ordering is specified by the
%                 algorithm input.
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
%    algorithm An optional parameter specifying the ordering into which the
%              rank is desired. Possible values are:
%              0 (The default if omitted or an empty matrix is passed) Use
%                colexicographic ordering via the combinatorial number
%                system that is described in Chapter 7.2.1.3 of [3].
%              1 Use lexicographic ordering, implementing the equation at
%                the top of page 181 of [1], which is essentially the
%                algorithm of [2]. If this algorithm is chosen, it is
%                required that the input n be provided.
%            n The number of items out of which the combinations are being
%              drawn. This input is only needed if algorithm=1.
%
%%OUTPUTS: theRank The rank of the combination in a lexicographic ordering
%                  of combinations, counting from 0.
%
%Colexicographic ordering has the advantage that when ranking a
%combination, one does not need to known n.
%
%EXAMPLE:
%We demonstrate that rankCombination is consistent with unrankCombination:
% n=10;
% m=7;
% startVal=0;
% firstElMostSig=false;
% algorithm=1;
% numCombos=binomial(n,m);
% for comboIdx=0:(numCombos-1)
%     curCombo=unrankCombination(comboIdx,n,m,firstElMostSig,startVal,algorithm);
%     theRank=rankCombination(curCombo,firstElMostSig,startVal,algorithm,n);
%     assert(theRank==comboIdx)
% end
%There should be no error (all assertions true), which indicates that the
%ranks are consistent with the unranked values.
%
%REFERENCES:
%[1] B. P. Buckles and M. Lybanon, "Algorithm 515: Generation of a vector
%    from the lexicographic index [G6]," ACM Transactions on Mathematical
%    Software, vol. 3, no. 2, pp. 180-182, Jun. 1977.
%[2] H. F. Walter, "Location of a vector in a lexicographically ordered
%    list," Communications of the ACM, vol. 6, no. 2, pg. 68, Feb. 1963.
%[3] D. E. Knuth, The Art of Computer Programming. Vol. 4, Fascicle 3:
%    Generating all Combinations and Partitions, Upper Saddle River, NJ:
%    Addison-Wesley, 2009.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(algorithm))
    algorithm=0;
end

if(nargin<3||isempty(startVal))
    startVal=0;
end

if(nargin<2||isempty(firstElMostSig))
    firstElMostSig=false;
end

switch(algorithm)
    case 0%Colexicographic order.
        theRank=rankColexCombination(combo,firstElMostSig,startVal);
    case 1%Lexicographic order.
        %we subtract 1, because the function ranks starting at 1.
        theRank=rankLexCombo(combo,n,firstElMostSig,startVal)-1;
    otherwise
        error('Unknown algorithm selected.')
end
end

function theRank=rankLexCombo(c,n,firstElMostSig,startVal)
%%RANKLEXCOMBINATION Obtain the lexicographic order (the rank) of a
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
%OUTPUTS: theRank The rank of the combination in a lexicographic ordering of
%                 combinations, counting from 1.
%
%This function implements the equation at the top of page 181 of [1]. It is
%essentially the same algorithm as in [2].
%
%REFERENCES:
%[1] B. P. Buckles and M. Lybanon, "Algorithm 515: Generation of a vector
%    from the lexicographic index [G6]," ACM Transactions on Mathematical
%    Software, vol. 3, no. 2, pp. 180-182, Jun. 1977.
%[2] H. F. Walter, "Location of a vector in a lexicographically ordered
%    list," Communications of the ACM, vol. 6, no. 2, pg. 68, Feb. 1963.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

p=length(c);

c=c-(startVal-1);
if(firstElMostSig)
    c=flipud(c(:));
end

x=[0;c(:)];

theRank=1;
for j=1:p
    %Note that x(j-1+1)+1 is not always less than or equal to x(j+1)-1, in
    %which case the for-loop is skipped.
    for k=(x(j-1+1)+1):1:(x(j+1)-1)
        theRank=theRank+binomial(n-k,p-j);
    end
end
end

function theRank=rankColexCombination(combo,firstElMostSig,startVal)
%%RANKCOLEXCOMBINATION Obtain the colexicographic order (the rank) of a
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
