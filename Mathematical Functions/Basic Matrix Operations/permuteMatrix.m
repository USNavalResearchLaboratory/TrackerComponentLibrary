function CPerm=permuteMatrix(C,order)
%%PERMUTEMATRIX Given an n1Xn2X...nS matrix, rearrange the dimensions of C
%      to be in the order specified by order. All elements of order must be
%      unique, real, positive, integer values from 1 to S. This functions
%      in the same manner as Matlab's permute function. The permutation
%      performed by this function is not done in-place. This type of
%      permutation is a type of tensor matrix transpose.
%
%INPUTS: C An n1Xn2X...XnS matrix.
%    order The length S order vector. This should hold integer values from
%          1 to S with no repeats.
%
%OUTPUTS: CPerm C with its dimensions permuted according to order.
%
%The algorithm is essentially brute force. It goes through each element in
%C and find the tuple of the associated element in CPerm to which the
%assignment must be performed.
%
%EXAMPLE:
% C=100*rand(3,4,5);
% order=[3,2,1];
% CPerm=permuteMatrix(C,order);
%One will ntoe that size(CPerm) is note [5,4,3].
%
%March 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

S=ndims(C);
numOrderEls=length(order);

if(isempty(order)&&~isempty(C))
   error('order cannot be empty with a non-empty C.') 
end

if(numOrderEls<S)
    error('The ordering vector must be >=ndims(C).')
end

if(any(order>numOrderEls)||any(order<1))
    error('')
end

if(isempty(C))
    CPerm=[];
    return;
end

nVals=size(C);

if(numOrderEls>S)
   nVals=[nVals,ones(1,numOrderEls-S)];
   S=numOrderEls;
end

totalEls=numel(C);

cumProdNew=zeros(S,1);
nValsNew=zeros(S,1);
invPerm=zeros(S,1);
newIdx=zeros(S,1);%Initialization
idx=zeros(S,1);%Initialization

iNew=order(1);
nValsNew(1)=nVals(iNew);

invPerm(iNew)=1;
cumProdNew(1)=1;

for i=2:S
    iNew=order(i);
    
    nValsNew(i)=nVals(iNew);
    invPerm(iNew)=i;

    cumProdNew(i)=cumProdNew(i-1)*nValsNew(i-1);
end

CPerm=zeros(nValsNew(:)');

newLinIdx=0;
linIdx=0;
while(1)
    %Assign the current tuple.
    CPerm(newLinIdx+1)=C(linIdx+1);
    
    %Move on to the next tuple.
    linIdx=linIdx+1;
    
    if(linIdx>=totalEls)
        break;
    end
    
    %The code below is similar to that in getNextTuple, because we are
    %updating the tuples for newIdx.
    curLevel=1;
    while(1)
        %Try incrementing the order at this level.
        idx(curLevel)=idx(curLevel)+1;
        if(idx(curLevel)<nVals(curLevel))
            iNew=invPerm(curLevel);
            newIdx(iNew)=newIdx(iNew)+1;
            newLinIdx=newLinIdx+cumProdNew(iNew);

            idx(1:(curLevel-1))=0;
            break;
        else
            %If the value is invalid, then just keep ascending and
            %adjust newLinIdx.
            iNew=invPerm(curLevel);

            newLinIdx=newLinIdx-newIdx(iNew)*cumProdNew(iNew);
            newIdx(iNew)=0;

            curLevel=curLevel+1;
            continue;
        end
    end
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