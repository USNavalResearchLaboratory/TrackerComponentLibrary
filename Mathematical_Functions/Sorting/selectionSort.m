function [A,idxList]=selectionSort(A,byCol,gtCompareFunc)
%%SELECTIONSORT Sort an array, a matrix, or a cell array in ascending order
%               using selection sort. This is a "straight" selection sort
%               algorithm, not a tree-based one. The algorthm is O(n^2)
%               best, average and worst case complexities, which is worse
%               than heapsort and quicksort. A custom comparison function
%               can be used so that one can, for example, sort strings or
%               other things. Changing the direction of the comparison
%               function (from greater than to less than) changes the
%               sorting order from ascending to descending. Note that
%               selectionSort is not a stable sorting algorithm, meaning
%               that the order of items having the same value might change.
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
%gtCompareFunc  A function handle that performs a greather-than comparison
%               of two entries in A. This lets one define custom comparison
%               operations. Providing a less-than comparison for this will
%               cause the list to be sorted in descending order rather than
%               ascending order. When sorting an array or matrix, the
%               the function handle takes inputs of the form
%               gtCompareFunc(A(:,i),A(:,j)) if byCol=true and with
%               reversed row and column indices if byCol=false. When
%               comparing cell arrays, it must handle inputs of the form
%               gtCompareFunc(A{i},A{j}). The default if this parameter is
%               omitted is @(x1,x2)(x1>x2);
%
%OUTPUTS: A The sorted array/ matrix/ cell array. If gtCompareFunc
%           performs a greater-than comparison, then it is in increasing
%           order. Otherwise it is in decreasing order.
%   idxList The indices of the original elements with respect to the sorted
%           order. For example, if the input A is an array, then
%           A(idxList) on the input A will give the sorted output A.
%
%The implementation of selectionSort is described in Chapter 5.2.3 of [1].
%
%Note that sorting large matrices can be slow as each element (so an entire
%row/ column) is copied during the search. Thus, if one is just sorting
%according to a particular row, it makes sense ot get the idxList for that
%one row and then use it to sort everything else.
%
%REFERENCES:
%[1] D. Knuth, The Art of Computer Programming: Sorting and Searching, 2nd
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
        
        numPoints=size(A,2);
    else
        numPoints=length(A);
    end
    
    %The index list is not used in the algorithm for sorting; it is just
    %computed in case the user wants it as a return variable.
    idxList=1:numPoints;

    if(~isACellArray)
    %If A is not a cell array.
        for j=numPoints:-1:2
            %Find the maximum value entry in the range 1:j.
            maxIdx=[];
            maxVal=-Inf;
            for i=j:-1:1
                if(~gtCompareFunc(maxVal,A(:,i)))
                    maxVal=A(:,i);
                    maxIdx=i;
                end
            end

            %Swap records idx and j.
            A(:,maxIdx)=A(:,j);
            A(:,j)=maxVal;
            
            %Record the swap in the index array
            temp=idxList(maxIdx);
            idxList(maxIdx)=idxList(j);
            idxList(j)=temp;
        end
        
        %Restore the orientation if it was supposed to be sorted by row.
        if(byCol==false)
            A=A';
        end
    else
    %We are here is A is a cell array. The only thing that changes is the
    %indexation of A.
        for j=numPoints:-1:2
            %Find the maximum value entry in the range 1:j.
            maxIdx=[];
            maxVal=A{j};
            for i=(j-1):-1:1
                if(~gtCompareFunc(maxVal,A{i}))
                    maxVal=A{i};
                    maxIdx=i;
                end
            end

            %Swap records idx and j.
            A{maxIdx}=A{j};
            A{j}=maxVal;
            
            %Record the swap in the index array
            temp=idxList(maxIdx);
            idxList(maxIdx)=idxList(j);
            idxList(j)=temp;
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
