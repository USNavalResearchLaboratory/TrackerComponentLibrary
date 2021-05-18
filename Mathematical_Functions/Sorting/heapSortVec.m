function [a,idxList]=heapSortVec(a,direction)
%%HEAPSORTVEC Sort a vector of scalar numeric values in ascending or
%          descending order using heap sort. Heap sort has a lower worst-
%          case performance bound O(n*log(n)) than quicksort O(n^2), but a
%          worse average case performance due to a multiplicative constant.
%          The best case complexity is also O(n*log(n)) Note that heapSort
%          is not a stable sorting algorithm, meaning that the order of
%          items having the same value might change. To sort more general
%          data types, consider the heapSort function.
%
%INPUTS: a The vector that is to be sorted.
% direction A value indicating the direction in which A should be sorted.
%          Possible values are boolean:
%          0 (The default if omitted or an empty matrix is passed) Sort in
%            ascending order.
%          1 Sort in descending order.
%
%OUTPUTS: a The input vector a with its elements sorted.
%   idxList The indices of the original vector with respect to the sorted
%           order.
%
%The algorithm for creating and updating the heap is generally based on the
%class implementation described in Chapter 6.4 of [1], though classes are
%not used here. The implementation using such a heap for sorting is
%described in Chapter 5.2.3 of [2].
%
%REFERENCES:
%[1] M.A.Weiss, Data Structures and Algorithm Analysis in C++, 2nd ed.
%    Reading, MA: Addison-Wesley, 1999.
%[2] D. Knuth, The Art of Computer Programming: Sorting and Searching, 2nd
%    ed. Reading, MA: Addison-Wesley, 1998, vol. 3.
%
%February 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %If an empty matrix is passed, return an empty matrix.
    if(isempty(a))
        idxList=[];
        return;
    end

    if(nargin<2||isempty(direction))
        direction=0;%Ascending order (the default).
    end
    
    %Step 1, build the heap.
    numInHeap=length(a);

    %The idx list plays no role in the algorithm and is just maintained
    %so that one can get a list of how the elements in A were
    %rearranged, if desired.
    idxList=1:numInHeap;

    idx=floor(numInHeap/2);
    while(idx>0)
        %Percolate the thing in index idx down into the heap.
        [a,idxList]=percolateDown(a,idxList,idx,numInHeap,direction);
        idx=idx-1; 
    end

    %Step 2, given the heap, the item at the beginning of the heap will be
    %the maximum item. Swap it with the item at the end, reduce the heap
    %size by 1 and reheapify thing. In the end, the array will be sorted in
    %ascending order.
    while(numInHeap>1)
        %Swap the first and last items.
        temp=a(1);
        a(1)=a(numInHeap);
        a(numInHeap)=temp;

        %Do the same for the index list.
        temp=idxList(1);
        idxList(1)=idxList(numInHeap);
        idxList(numInHeap)=temp;

        numInHeap=numInHeap-1;
        %Reheapify
        [a,idxList]=percolateDown(a,idxList,1,numInHeap,direction);
    end
end

function [A,idxList]=percolateDown(A,idxList,idx,numInHeap,direction)
    %Note the idxList plays no role in the algorithm; it is just used to
    %keep track of how the indices changed, in case the user wants that
    %information.

    temp=A(idx);
    tempIdx=idxList(idx);
    if(direction==0)%Ascending direction    
        while(2*idx<=numInHeap)
            child=2*idx;

            if(child~=numInHeap&&A(child+1)>A(child))
               child=child+1; 
            end

            if(A(child)>temp)
                A(idx)=A(child);
                idxList(idx)=idxList(child);
            else
                break;
            end

            idx=child; 
        end
    else%Descending direction
         while(2*idx<=numInHeap)
            child=2*idx;

            if(child~=numInHeap&&A(child+1)<A(child))
               child=child+1; 
            end

            if(A(child)<temp)
                A(idx)=A(child);
                idxList(idx)=idxList(child);
            else
                break;
            end

            idx=child; 
         end
    end

    A(idx)=temp;
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
