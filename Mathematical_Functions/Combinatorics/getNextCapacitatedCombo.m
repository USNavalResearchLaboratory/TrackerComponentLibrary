function [theCombo,recurData]=getNextCapacitatedCombo(param1,numBins,totalItems)
%%GETNEXTCAPACITATEDCOMBO Get the next combination of placing totalItems
%               items into numBins unlabeled bins where each bin can hold
%               numPerBin items. If numPerBin=1, then this is just
%               generating standard combinations of numBins items taken
%               from a total of totalItems. The capacitated combinations
%               are equivalent to putting totalItems labeled balls into
%               numBins unlabeled bins such that each bin holds numPerBin
%               items and the ordering of the items in the bins does not
%               matter.
%
%INPUTS: param1 If the first combination is desired, then more than one
%               input must be given and param1=numPerBin, the number of
%               items in each bin. Otherwise, param1=recurData, the second
%               output of the function from the last call.
%       numBins This input should only be provided if the first capacitated
%               combination in the sequence is desired. This is the number
%               of bins (of capacity numPerBin) into which things can be
%               placed.
%    totalItems This input can only be given when the first combination is
%               desired. This is the total number of items that can be
%               assigned. If numBins is provided, but this is omitted or an
%               empty matrix is passed, then a value of numBins*numPerBin
%               is used, the minimum allowed value.
%
%OUTPUTS: theCombo The numPerBinXnumBins next capacitated combination. If
%                  the final capacitated combination has been passed, then
%                  an empty matrix is returned.
%        recurData A data structure that can be passed back to this
%                  function to get the next combination.
%
%The algorithm is the same as that used in genAllCapacitatedCombos. See the
%comments to genAllCapacitatedCombos for more information. The total number
%of capaciatted combinations is given by numCapacitatedCombos.
%
%EXAMPLE:
%This example just shows that the capacitated combinations produced by this
%function come in the same order as those produced by the function
%genAllCapacitatedCombos. When run, it collects the combinations by
%repeatedly calling this function and then asserts that they are all equal
%what genAllCapacitatedCombos returns. The result shoud be true, so there
%should be no error.
% numPerBin=2;
% numBins=4;
% totalItems=4*2+2;
% 
% totalCapCombos=numCapacitatedCombos(numPerBin,numBins,totalItems);
% allCapCombos=zeros(numPerBin,numBins,totalCapCombos);
% 
% [theCombo,recurData]=getNextCapacitatedCombo(numPerBin,numBins,totalItems);
% allCapCombos(:,:,1)=theCombo;
% 
% for comboIdx=2:totalCapCombos
%     [theCombo,recurData]=getNextCapacitatedCombo(recurData);
%     allCapCombos(:,:,comboIdx)=theCombo;
% end
% 
% allBatchCombos=genAllCapacitatedCombos(numPerBin,numBins,totalItems);
% assert(all(allCapCombos(:)==allBatchCombos(:)))
%The assertion should not raise an error.
%
%July 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin>=2)
    if(nargin<3||isempty(totalItems))
        totalItems=param1*numBins;
    end

    %If this is the first capacitated combination, then create the
    %recurData structure for the first time.
    numPerBin=param1;
    
    %The total number of items to assign in combinations across the bins.
    n=numBins*numPerBin;

    %Save the values in recurData that do not change with the iterations.
    recurData=[];
    recurData.numPerBin=numPerBin;
    recurData.numBins=numBins;
    recurData.n=n;
    recurData.totalItems=totalItems;
    
    %Special case for standard combinations.
    if(numPerBin==1)
        curIdxCombos=1:numBins;
        recurData.curIdxCombos=curIdxCombos;
        theCombo=curIdxCombos;
        recurData.curBin=1;
        recurData.curSet=1;
        recurData.curCombos=[];
        recurData.idxLists=[];
        recurData.backtracking=[];
        
        I=1:n;
        recurData.I=getNextCombo(I,totalItems,1);
        return
    end
    
    curIdxCombos=zeros(numPerBin,numBins);
    curCombos=zeros(numPerBin-1,numBins);

    %This holds the indices that can be assigned in each bin and in subsequent
    %bins.
    idxLists=zeros(n,numBins);
    idxLists(:,1)=1:n;

    curSet=1;
    curBin=1;
    backtracking=false;
    %When n~=totalItems, after each base capacitated combination is
    %generated, we have to go through all combinations of indices that make
    %up the combinations. This is the current index combination. When
    %empty, it means that we haven't gotten a base capacitated combination
    %yet.
    I=[];
    recurData.I=I;
else
    recurData=param1;
    
    if(isempty(recurData))
       theCombo=[];
       return;
    end

    %Extract the values that change on each iteration.
    curIdxCombos=recurData.curIdxCombos;
    curCombos=recurData.curCombos;
    idxLists=recurData.idxLists;
    curSet=recurData.curSet;
    curBin=recurData.curBin;
    backtracking=recurData.backtracking;
    I=recurData.I;
end

%Get the values that do not change between calls.
numPerBin=recurData.numPerBin;
numBins=recurData.numBins;
n=recurData.n;
totalItems=recurData.totalItems;
numPerBinM1=numPerBin-1;

if(~isempty(I))
    %Just get the next combination of indices to the capcitated combination
    %and return that. 
    theCombo=reshape(I(curIdxCombos),[numPerBin,numBins]);

    %If past the final combination, then I is an empty matrix and it will
    %know to go onto the next capacitated base combination.
    I=getNextCombo(I,totalItems,1);
    recurData.I=I;
    return
elseif(numPerBin==1)
    %The special case of standard combos.
    theCombo=[];
    return;
end

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
            %If we are in the final bin, then save the combo and move to
            %the previous bin. In practice, subsequent calls to this
            %function will nto go to the previous bin until all
            %combinations of index assignments (via I) are assigned.
            I=1:n;
            theCombo=curIdxCombos;
            curSet=curSet+1;
            
            backtracking=true;
            curBin=curBin-1;
            
            %Save the values that change during iteration to recurData
            %before returning.
            recurData.curIdxCombos=curIdxCombos;
            recurData.curCombos=curCombos;
            recurData.idxLists=idxLists;
            recurData.curSet=curSet;
            recurData.curBin=curBin;
            recurData.backtracking=backtracking;
            recurData.I=getNextCombo(I,totalItems,1);
            return;
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

%If we get here, we passed the final combiantion.
theCombo=[];
recurData=[];
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
