function tuples=genAllUnorderedTuples(numEls,maxVal,firstMostSig)
%%GENALLUNORDEREDTUPLES  Unordered tuples are tuples where the ordering of
%           the elements in the tuple does not matter. Thus, the tuples
%           1,1,2,3 and 3,1,2,1 represent the same value. This function
%           generates all length-numEls tuples with element values ranging
%           from 0 to maxVal.
%
%INPUTS: numEls The number of elements in each tuple.
%       maxVal The maximum value that each digit of the tuple can take
%              (>=0). This corresponds to the base of each digit -1.
%  firstIsMostSig This is a boolean variable indicating whether tuple(1) is
%              the most significant digit (or whether tuple(N) is the most
%              significant). This affects the ordering of the sequence that
%              produces the unordered tuples. The default if omitted or an
%              empty matrix is passed is false.
%
%OUTPUTS: tuples AnumElsXnumTUples matrix of all of the tuples.
%
%This function kist calls getNextUnorderedTuple in a loop. The total number
%of unordered tuples is agiven by the numUnorderedTuples function.
%
%EXAMPLE:
% tuples=genAllUnorderedTuples(3,4)
%One will get
%tuple=[0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 2;
%       0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 2, 2, 2;
%       0, 0, 1, 1, 1, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2;
%       0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(firstMostSig))
    firstMostSig=false;
end

numTup=numUnorderedTuples(numEls,maxVal);
tuples=zeros(numEls,numTup);
tuples(:,1)=zeros(numEls,1);
for curTup=2:numTup
    tuples(:,curTup)=getNextUnorderedTuple(tuples(:,curTup-1),maxVal,firstMostSig);
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
