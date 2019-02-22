function x=runLenDecode(vals,numReps)
%%RUNLENDECODE Given a run-length encoded message --that is a set of values
%              and a vector saying how many times each value is repeated,
%              obtain the original string. That is, decode the run-length
%              encoded string.
%
%INPUTS: vals A set of values.
%     numReps The number of times each element in vals is repeated to
%             obtain the original string x. These values must be >=0.
%
%OUTPUTS: x The reconstructed vector as a column vector. The number of
%           elements in the vector equals sum(numReps).
%
%Run length encoding can be a simple, useful data compression method when
%dealing with data containing many consecutive repeated values, such as a
%black and white bitmap image.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numUniqueChars=length(vals);
totalChars=sum(numReps);

x=zeros(totalChars,1);
numInX=0;
for i=1:numUniqueChars
    numRepeats=numReps(i);
    
    x((numInX+1):(numInX+numRepeats))=vals(i);
    numInX=numInX+numRepeats;
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
