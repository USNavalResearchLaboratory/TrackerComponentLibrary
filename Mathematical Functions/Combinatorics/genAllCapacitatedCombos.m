function allCombos=genAllCapacitatedCombos(numPerBin,numBins,totalItems)
%%GENALLCAPACITATEDCOMBOS Given totalItems items, generate all combinations
%               combinations of placing the items into numBins unlabeled
%               bins where each bin can hold numPerBin items. If
%               numPerBin=1, then this is just generating standard
%               combinations of numBins items taken from a total of
%               totalItems. The capacitated combinations are equivalent to
%               putting totalItems labeled balls into numBins unlabeled
%               bins such that each bin holds numPerBin items and the
%               ordering of the items in the bins does not matter.
%
%INPUTS: numPerBin The scalar integer capacity of the bins; numPerBin>=1.
%          numBins The scalar integer number of bins; numBins>=1.
%       totalItems The scalar integer total number of items to place in the
%                  bins; totalitems>=numBins*numPerBin. If omitted or an
%                  empty matrix is passed, then the default of 
%                  totalItems=numBins*numPerBin is used.
%
%OUTPUTS: allCombos The numPerBinXnumBinsXnumCapCombos set of capacitated
%                   combinations. allCombos(:,:,i) is the ith capacitated
%                   combination. The rows are the items in each bin and the
%                   columns are the bins. Items are numbered from 1 to
%                   totalItems.
%
%If we assume that totalItems=n=numPerBin*numBins, then we can generate the
%combinations by fixing the first element of each bin to the smallest index
%unassigned to any previous bin in the current combination being
%constructed. Thus, the first bin will always contain item 1. The second
%bin will always contain the smallest unassigned item that is not in the
%first bin and so on. Thus, we have to go through standard
%(n-(k-1)*numPerBin-1) choose (numPerBin-1) combinations of the remaining
%indices for the kth bin conditioned on the previous bins. This function
%implements a non-recursive algorithm for generating such combinations.
%
%In the general case, if totalItems>n, then we have to additionally go
%through all possible combinations of n items chosen from totalItems. This
%means generating all combinations as if totalItems=n and then remapping
%the indices for each subset of n items from totalItems.
%
%A total of numCapacitatedCombos(numPerBin,numBins,totalItems) capacitated
%combinations will be produced. Note that the number of capacitated
%combinations can grow rapidly and one might prefer visiting them
%sequentially using getNextCapacitatedCombo rather than generating all of
%them at once.
%
%EXAMPLE:
%Here we put six items into two bins of capacity 3:
% allCombos=genAllCapacitatedCombos(3,2,6)
%One wil get 10 combinations:
% allCombos(:,:,1)=[1, 4;
%                   2, 5;
%                   3, 6];
% allCombos(:,:,2)=[1, 3;
%                   2, 5;
%                   4, 6];
% allCombos(:,:,3)=[1, 3;
%                   2, 4;
%                   5, 6];
% allCombos(:,:,4)=[1, 3;
%                   2, 4;
%                   6, 5];
% allCombos(:,:,5)=[1, 2;
%                   3, 5;
%                   4, 6];
% allCombos(:,:,6)=[1, 2;
%                   3, 4;
%                   5, 6];
% allCombos(:,:,7)=[1, 2;
%                   3, 4;
%                   6, 5];
% allCombos(:,:,8)=[1, 2;
%                   4, 3;
%                   5, 6];
% allCombos(:,:,9)=[1, 2;
%                   4, 3;
%                   6, 5];
% allCombos(:,:,10)=[1, 2;
%                    5, 3;
%                    6, 4];
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The total number of items to assign in combinations across the bins.
n=numBins*numPerBin;

if(nargin<3||isempty(totalItems))
    totalItems=n;
end

if(n>totalItems)
    error('totalItems cannot be less than numBins*numPerBin')
end

%Special case for standard combinations.
if(numPerBin==1)
    allCombos=genAllCombinations(totalItems,numBins,1,1);
    numCopies=size(allCombos,2);
    allCombos=reshape(allCombos,[numPerBin,numBins,numCopies]);
    return
end

%The total number of items to assign in combinations.
n=numBins*numPerBin;

%The total number of capacitated combinations if n=totalItems.
numCombosNoExtra=numCapacitatedCombos(numPerBin,numBins);

%If n~=totalItems, then this is the number of variants of each capacitated
%combination that will have to be generated.
numCopies=binomial(totalItems,n);

%To store the results.
allCombos=zeros(numPerBin,numBins,numCopies,numCombosNoExtra);

numPerBinM1=numPerBin-1;

curIdxCombos=zeros(numPerBin,numBins);
curCombos=zeros(numPerBin-1,numBins);

%This holds the indices that can be assigned in each bin and in subsequent
%bins.
idxLists=zeros(n,numBins);
idxLists(:,1)=1:n;

curSet=1;
curBin=1;
backtracking=false;
while(curBin>0)
    %This is the number of elements that can be assigned in the current
    %bin including making a hard assignment on the first index.
    numInIdx=n-(curBin-1)*numPerBin;

    if(backtracking==false)
        %We are just entering this bin.
        %The first element in each bin is fixed to the first index in the
        %idxList. 
        curIdxCombos(1,curBin)=idxLists(1,curBin);

        %In this level, we have to go through all combinations of the
        %remaining n-curBin*(numPerBin-1); unassigned elements of
        %idxLists(:,curBin) in this bin. Each time an assignment is made in
        %this bin, those elements are removed from idxLists going into the
        %next higher level so that we never get repeats.
        curIdxCombos(2:numPerBin,curBin)=idxLists(2:numPerBin,curBin);
        %Current combination of the indices in idxLists in this bin.
        curCombos(:,curBin)=1:numPerBinM1;

        if(curBin==numBins)
            %If we are in the final bin, then save the combo with all
            %variants of indices and move to the previous bin.
            
            I=1:n;
            for curCopy=1:numCopies
                allCombos(:,:,curCopy,curSet)=I(curIdxCombos);

                I=getNextCombo(I,totalItems,1);
            end

            curSet=curSet+1;
            
            backtracking=true;
            curBin=curBin-1;
            continue;
        else
            %Otherwise, set idxLists for the next level and move to the
            %next bin.numInIdxNext
            numInIdxNext=numInIdx-numPerBin;
            idxLists(1:numInIdxNext,curBin+1)=idxLists((numPerBin+1):numInIdx,curBin);

            curBin=curBin+1;
        end
    else
        %If here, we are backtracking.

        %We subtract 1 because we have made a hard assignment on the first
        %index.
        nextCombo=getNextCombo(curCombos(:,curBin),numInIdx-1,1);
        if(isempty(nextCombo))
            %Continue backtracking.
            curBin=curBin-1;
            continue;
        end
        %Otherwise, record the current combination.
        curCombos(:,curBin)=nextCombo;
        
        %Record the indices for the current combination.
        nextCombo=nextCombo+1;
        curIdxCombos(2:numPerBin,curBin)=idxLists(nextCombo,curBin);
        
        %Now, set the values in idxLists at the next bin level to those
        %that are in idxLists but not in nextCombo. We note that the
        %indices are given in increasing order and that the first index is
        %skipped because we know it is assigned.
        curAssignedIdx=2;
        curNextIdx=1;
        for curListIdx=2:numInIdx
            if(idxLists(curListIdx,curBin)==curIdxCombos(curAssignedIdx,curBin))
                curAssignedIdx=curAssignedIdx+1;
                if(curAssignedIdx>numPerBin)
                    %Then all of the rest of the indices are assigned.
                    numInIdxNext=numInIdx-numPerBin;
                    idxLists(curNextIdx:numInIdxNext,curBin+1)=idxLists((curListIdx+1):numInIdx,curBin);
                    break;
                end
            else
                idxLists(curNextIdx,curBin+1)=idxLists(curListIdx,curBin);
                curNextIdx=curNextIdx+1;
            end
        end

        %Go to the next bin.
        curBin=curBin+1;
        backtracking=false;
    end
end

allCombos=allCombos(:,:,:);

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
