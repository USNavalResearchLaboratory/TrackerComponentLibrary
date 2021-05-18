function numTup=numUnorderedTuples(numEls,maxVal)
%%NUMUNORDEREDTUPLES Get the total number of unordered tuples. Unordered
%           tuples are tuples where the ordering of the elements in the
%           tuple does not matter. Thus, the tuples 1,1,2,3 and 3,1,2,1
%           represent the same value.
%
%INPUTS: numEls The number of elements in the tuple. numEls>=1.
%        maxVal The maximum value that each element of the tuple can take.
%               It goes from 0 to this value, so maxVal+1 is the base used
%               for counting each element.
%
%OUTPUTS: numTup The number of unordered tuples.
%
%This is a counting problem. The solution is expressed in terms of a
%binomial value. Unordered tuples are related to multiset combinations.
%Indeed, this function produces the same result as
% m=ones(maxVal+1,1)*numEls;
% numMultisetCombos(m,numEls)
%but is simpler.
%
%EXAMPLE:
%Here, we generate all unordered tuples and verify that the number agrees
%with that is returned by this function
% numEls=5;
% maxVal=3;
% tuple=zeros(numEls,1);
% count=0;
% while(~isempty(tuple))
% 	count=count+1;
%     tuple=getNextUnorderedTuple(tuple,maxVal,true);
% end
% count
% numTup=numUnorderedTuples(numEls,maxVal)
%
%November 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numTup=binomial(numEls+maxVal,maxVal);

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
