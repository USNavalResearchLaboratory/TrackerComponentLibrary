function [vals,numReps]=runLenEncode(x)
%%RUNLENENCODE Given a vector of values, perform run length encoding. This
%              means returning a vector of values with no consecutive
%              repeats and providing a vector of how many consecutive
%              repeats each value has.
%
%INPUTS: x A numValsX1 or 1XnumVals vector.
%
%OUTPUTS: vals The values in x with no consecutive repeated elements.
%      numReps The number of times each element in vals is repeated to
%              obtain the original string x.
%
%Run length encoding can be a simple, useful data compression method when
%dealing with data containing many consecutive repeated values, such as a
%black and white bitmap image.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Allocate the maximum possible amount of space for the variables.
numVals=length(x);

%Allocate
vals=zeros(numVals,1);
numReps=zeros(numVals,1);

curIdx=1;
vals(curIdx)=x(1);
numReps(curIdx)=1;
for i=2:numVals
    if(x(i-1)==x(i))
        numReps(curIdx)=numReps(curIdx)+1;
    else
        curIdx=curIdx+1;
        vals(curIdx)=x(i);
        numReps(curIdx)=1;
    end
end

%Shrink to fit the actual number of characters present.
vals=vals(1:curIdx);
numReps=numReps(1:curIdx);

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
