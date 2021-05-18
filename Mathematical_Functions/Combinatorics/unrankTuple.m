function tuple=unrankTuple(rank,maxVals,areProdVals)
%%UNRANKTUPLE Obtain the tuple corresponding to its order in the sequence
%             of tuples given a certain set of maximum values for each
%             digit. These are generated with the first digit being the
%             LEAST significant.
%
%INPUTS: rank The rank of the tuple in the order of increasing tuples where
%             the first digit is the most significant. This ranges from 0
%             to the total possible number of tuples-1. The total number of
%             tuples depends on maxVals.
%      maxVals An NX1 or 1XN vector of the maximum value that each digit
%             of the tuple can take. Alternatively, if many calls are to be
%             made to this function, the function will be faster is one
%             sets the input areProdVals=true in which case, this input
%             should be [1;cumprod(maxVals+1)] instead of maxVals.
% areProdVals An optional boolean parameter specifying whether maxVals are
%             the maximum values for each digit or whether they are product
%             values as described above. The default if omitted or an empty
%             matrix is passed is false.
%
%OUTPUTS: tuple The NX1 tuple corresponding to the given rank.
%
%The concept is similar to how the index2NDim function works. That function
%unranks tuples, but it assumes that the digits start from 1 and not zero
%and that the rank starts from 1, not 0.
%
%EXAMPLE:
%Here, we show how to generate the same output using the two different
%input formulations of this function.
% maxVals=[6;3;4;2];
% idx=[149;18;9;3;0];
% areProdVals=false;
% indices0=unrankTuple(idx,maxVals,areProdVals)
% 
% maxVals=[1;cumprod(maxVals+1)];
% areProdVals=true;
% indices1=unrankTuple(idx,maxVals,areProdVals)
%One will see that indices0=indices1 and that
%indices0=[2, 4, 2, 3, 0;
%          1, 2, 1, 0, 0;
%          0, 0, 0, 0, 0;
%          1, 0, 0, 0, 0]
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(areProdVals))
    areProdVals=false;
end

numIdx=length(rank);

%The indices basically form a counting system where the bases of each
%position are determined by dim. Matlab makes the first dimension the
%least significant. Thus, finding the indices is done in the same manner as
%one would decompose a number into digits or bits.
if(areProdVals==true)
    numDim=length(maxVals)-1;
    maxVal=maxVals(end);
else
    numDim=length(maxVals);
    maxVals=[1;cumprod(maxVals+1)];
    maxVal=maxVals(numDim+1);
end

%If the index is above the maximum number of items.
if(any(rank>=maxVal))
   tuple=[];
   return;
end

tuple=zeros(numDim,numIdx);
rank=rank(:)';%Make it a row vector
for curIndic=numDim:-1:1
    %The number of complete multiples
    wholeVal=floor(rank/maxVals(curIndic));
    tuple(curIndic,:)=wholeVal;
    rank=rank-wholeVal*maxVals(curIndic);
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
