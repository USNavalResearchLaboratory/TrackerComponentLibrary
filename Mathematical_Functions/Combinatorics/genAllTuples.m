function theTuples=genAllTuples(maxVals,firstIsMostSig)
%%GENALLTUPLES Generate all tuples in a counting series. Each index of the
%              tuple counts from 0 to the corresponding value in maxVals.
%              This is just a method of counting uses different bases for
%              each digit.
%
%INPUTS: maxVals An NX1 or 1XN vector of the bases of each of the N digits
%                in the tuples minus 1. All elements must be >=0.
% firstIsMostSig This is a boolean variable indicating whether tuple(1) is
%                the most significant digit (or whether tuple(N) is the
%                most significant). The default if this parameter is
%                omitted or an empty matrix is passed is true.
%
%OUTPUTS: theTuples An NXnumTuples(maxVals) matrix of all of the tuples.
%                The ordering is such that theTuples(1,k) is the most
%                significant digit of the kth tuple and theTuples(N,k) is
%                the least significant digit of the kth tuple.
%
%This function just calls getNextTuple in a loop. The numTuples function
%returns the total number of tuples.
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(firstIsMostSig))
    firstIsMostSig=true; 
end

N=length(maxVals);
tupleCount=numTuples(maxVals);
theTuples=zeros(N,tupleCount);
theTuples(:,1)=getNextTuple(N);
for k=2:tupleCount
    theTuples(:,k)=getNextTuple(theTuples(:,k-1),maxVals,firstIsMostSig);
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
