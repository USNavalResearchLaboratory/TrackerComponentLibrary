function [A,idxList]=heapSort(A,byCol,gtCompareFunc)
%%HEAPSORT Sort an array, a matrix, or a cell array in ascending order
%          using a heap sort algorithm. A custom comparison function can
%          be used so that one can, for example, sort strings or other
%          things. Changing the direction of the comparison function (from
%          greater than to less than) changes the sorting order from
%          ascending to descending. Heap sort has a lower worst-case
%          performance bound O(n*log(n)) than quicksort O(n^2), but a
%          worse average case performance due to a multiplicative
%          constant. The best case complexity is also O(n*log(n)) Note
%          that heapSort is not a stable sorting algorithm, meaning that
%          the order of items having the same value might change. When only
%          sorting a vector of scalar values, consider the heapSortVec
%          function.
%
%INPUTS: A  An array, 2D matrix, or linear cell array that is to be sorted.
%           When given a matrix, the comparison function gtCompareFunc must
%           be provided so that it is clear how the columns are compared
%           (or the rows if byCol is false).
%     byCol A boolean value indicating whether sorting should be performed
%           by row or by column. This must be specified when a matrix is
%           passed as as a linear array could be confused with a matrix
%           with just one row/column. The default if this parameter is
%           omitted and the input is 1-dimensional is whatever would sort
%           over the 1D array. If the input is 2D, then the default is
%           true. This parameter is not used when cell arrays are passed
%           (an empty matrix can be passed for this) as only linear cell
%           arrays are supported, so no ambiguity between 1D and 2D inputs
%           would exist.
% gtCompareFunc A function handle that performs a greater-than comparison
%           of two entries in A. This lets one define custom comparison
%           operations. Providing a less-than comparison for this will
%           cause the list to be sorted in descending order rather than
%           ascending order. When sorting an array or matrix, the
%           function handle takes inputs of the form
%           gtCompareFunc(A(:,i),A(:,j)) if byCol=true and with reversed
%           row and column indices if byCol=false. When comparing cell
%           arrays, it must handle inputs of the form
%           gtCompareFunc(A{i},A{j}). The default if this parameter is
%           omitted is @(x1,x2)(x1>x2);
%
%OUTPUTS: A The sorted array/ matrix/ cell array. If gtCompareFunc
%           performs a greater-than comparison, then it is in increasing
%           order. Otherwise it is in decreasing order.
%   idxList The indices of the original elements with respect to the sorted
%           order. For example, if the input A is an array, then
%           A(idxList) on the input A will give the sorted output A.
%
%The algorithm for creating and updating the heap is generally based on the
%class implementation described in Chapter 6.4 of [1], though classes are
%not used here. The implementation using such a heap for sorting is
%described in Chapter 5.2.3 of [2].
%
%Note that sorting large matrices can be slow as each element (so an entire
%row/ column) is copied during the search. Thus, if one is just sorting
%according to a particular row, it makes sense ot get the idxList for that
%one row and then use it to sort everything else.
%
%REFERENCES:
%[1] M.A.Weiss, Data Structures and Algorithm Analysis in C++, 2nd ed.
%    Reading, MA: Addison-Wesley, 1999.
%[2] D. Knuth, The Art of Computer Programming: Sorting and Searching, 2nd
%    ed. Reading, MA: Addison-Wesley, 1998, vol. 3.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %If an empty matrix is passed, return an empty matrix.
    if(isempty(A))
        idxList=[];
        return;
    end

    if(nargin<2||isempty(byCol))
        %If the input is 1-dimensional, then make the byCol parameter go
        %over whichever dimension would sort it. Otherwise, make it go over
        %columns.
        if(size(A,1)==1)
            byCol=true;
        elseif(size(A,2)==1)
            byCol=false;
        else
            byCol=true;
        end
    end

    if(nargin<3||isempty(gtCompareFunc))
        gtCompareFunc=@(x1,x2)(x1>x2);
    end
    
    isACellArray=isa(A,'cell');
    
    if(~isACellArray)
        %If A is supposed to be sorted by row.
        if(byCol==false)
            A=A';
        end
       
        %Step 1, build the heap.
        numInHeap=size(A,2);
        
        %The idx list plays no role in the algorithm and is just maintained
        %so that one can get a list of how the elements in A were
        %rearranged, if desired.
        idxList=1:numInHeap;
        
        idx=floor(numInHeap/2);
        while(idx>0)
            %Percolate the thing in index idx down into the heap.
            [A,idxList]=percolateDown(A,idxList,idx,gtCompareFunc,numInHeap);
            idx=idx-1; 
        end

        %Step 2, given the heap, the item at the beginning of the heap will be
        %the maximum item. Swap it with the item at the end, reduce the heap
        %size by 1 and reheapify thing. In the end, the array will be sorted in
        %ascending order.
        while(numInHeap>1)
            %Swap the first and last items.
            temp=A(:,1);
            A(:,1)=A(:,numInHeap);
            A(:,numInHeap)=temp;
            
            %Do the same for the index list.
            temp=idxList(1);
            idxList(1)=idxList(numInHeap);
            idxList(numInHeap)=temp;
            
            numInHeap=numInHeap-1;
            %Reheapify
            [A,idxList]=percolateDown(A,idxList,1,gtCompareFunc,numInHeap);
        end

        %Restore the orientation if it was supposed to be sorted by row.
        if(byCol==false)
            A=A';
        end
    else
        %If we are here, then A is a cell array. Everything is pretty much
        %the same as above, but the indexation in Matlab has to change.
        
        %Step 1, build the heap.
        numInHeap=length(A);
        idxList=1:numInHeap;
        idx=floor(numInHeap/2);
        while(idx>0)
            %Percolate the thing in index idx down into the heap.
            [A,idxList]=percolateDownCell(A,idxList,idx,gtCompareFunc,numInHeap);
            idx=idx-1; 
        end

        %Step 2, given the heap, extract the elements in order.
        while(numInHeap>1)
            %Swap the first and last items.
            temp=A{1};
            A{1}=A{numInHeap};
            A{numInHeap}=temp;
            
            %Do the same for the index list.
            temp=idxList(1);
            idxList(1)=idxList(numInHeap);
            idxList(numInHeap)=temp;
            
            numInHeap=numInHeap-1;

            %Reheapify
            [A,idxList]=percolateDownCell(A,idxList,1,gtCompareFunc,numInHeap);
        end
    end
end

function [A,idxList]=percolateDown(A,idxList,idx,gtCompareFunc,numInHeap)
    %Note the idxList plays no role in the algorithm; it is just used to
    %keep track of how the indices changed, in case the user wants that
    %information.

    temp=A(:,idx);
    tempIdx=idxList(idx);
    while(2*idx<=numInHeap)
        child=2*idx;

        if(child~=numInHeap&&gtCompareFunc(A(:,child+1),A(:,child)))
            child=child+1; 
        end

        if(gtCompareFunc(A(:,child),temp))
            A(:,idx)=A(:,child);
            idxList(idx)=idxList(child);
        else
            break;
        end

        idx=child; 
    end
    A(:,idx)=temp;
    idxList(idx)=tempIdx;
end

function [A,idxList]=percolateDownCell(A,idxList,idx,gtCompareFunc,numInHeap)
%This is the same as percolateDown, but is for cell arrays.
    temp=A{idx};
    tempIdx=idxList(idx);
    while(2*idx<=numInHeap)
        child=2*idx;

        if(child~=numInHeap&&gtCompareFunc(A{child+1},A{child}))
            child=child+1; 
        end

        if(gtCompareFunc(A{child},temp))
            A{idx}=A{child};
            idxList(idx)=idxList(child);
        else
            break;
        end

        idx=child; 
    end
    A{idx}=temp;
    idxList(idx)=tempIdx;
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
