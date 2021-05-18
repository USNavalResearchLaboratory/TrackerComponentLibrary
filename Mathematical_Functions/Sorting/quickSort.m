function [A,idxList]=quickSort(A,byCol,gtCompareFunc)
%%QUICKSORT Sort an array, a matrix, or a cell array in ascending order
%           using a non-recursive implementation of the quicksort
%           algorithm. The quicksort algorithm tends to have better average
%           performance than other techniques, being O(n*log(n)) averge and
%           best case, with a small multiplicative constant, but has O(n^2)
%           worst case performance, which is worse than heapsort.  A custom
%           comparison function can be used so that one can, for example,
%           sort strings or other things. Changing the direction of the
%           comparison function (from greater than to less than) changes
%           the sorting order from ascending to descending. Note that
%           quickSort is not a stable sorting algorithm, meaning that the
%           order of items having the same value might change.
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
%The quicksort algorithm is based on the description given in Chapter 7 of
%[1]. However, it has been modified to eliminate the recursion.
%
%Note that sorting large matrices can be slow as each element (so an entire
%row/ column) is copied during the search. Thus, if one is just sorting
%according to a particular row, it makes sense ot get the idxList for that
%one row and then use it to sort everything else.
%
%REFERENCES:
%[1] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein,
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
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
    
    %stack(:,i) holds the lower and upper bounds of the partiion being
    %sorted at the ith level. The upper bound on the stack space in
    %the worst-case scenario is where every single point visited
    %creates two partitions. However, since each time it splits, it
    %deletes itself, the limit on the memory is ceil(numPoints/2) and
    %not numPoints.
    stack=zeros(2,ceil(numPoints/2));

    p=1;
    r=numPoints;
    stackIdx=1;
    stack(:,stackIdx)=[p;r];
    while(stackIdx>0)
        %Pop the p and r values off of the stack.
        p=stack(1,stackIdx);
        r=stack(2,stackIdx);
        stackIdx=stackIdx-1;

        %Partition the array, getting the correct location of the pivor
        %element.
        if(~isACellArray)
            [q,A,idxList]=partition(A,idxList,p,r,gtCompareFunc);
        else
            [q,A,idxList]=partitionCell(A,idxList,p,r,gtCompareFunc);
        end

        %If there are elements on the left side, then push them onto
        %the stack.
        if(p<q-1)
            stackIdx=stackIdx+1;
            stack(:,stackIdx)=[p;q-1];
        end

        %If there are elements on the right side, then push them onto
        %the stack.
        if(q+1<r)
            stackIdx=stackIdx+1;
            stack(:,stackIdx)=[q+1;r];
        end
    end
    
    %Restore the orientation if it was supposed to be sorted by row.
    if(~isACellArray&&byCol==false)
        A=A';
    end
end

function [q,A,idxList]=partition(A,idxList,p,r,gtCompareFunc)
%Rearrange the subarray A(p:r) in place.
    x=A(:,r);
    i=p-1;
    for j=p:(r-1)
        if(~gtCompareFunc(A(:,j),x))
            i=i+1; 
            %Swap A(:,i) and A(:,j)
            temp=A(:,i);
            A(:,i)=A(:,j);
            A(:,j)=temp;
            
            %Swap the indices in idxList to keep track of changes in the
            %ordering in A.
            temp=idxList(i);
            idxList(i)=idxList(j);
            idxList(j)=temp;
        end
    end

    %Swap A(:,i+1) and A(:,r)
    temp=A(:,i+1);
    A(:,i+1)=A(:,r);
    A(:,r)=temp;
    
    %Record the swap in the index list.
    temp=idxList(i+1);
    idxList(i+1)=idxList(r);
    idxList(r)=temp;
    
    %The new pivot point to return.
    q=i+1;
end

function [q,A,idxList]=partitionCell(A,idxList,p,r,gtCompareFunc)
%This is the same as the partition function, but has been modified to index
%cell arrays.
    x=A{r};
    i=p-1;
    for j=p:(r-1)
        if(~gtCompareFunc(A{j},x))
            i=i+1; 
            %Swap A{i} and A{j}
            temp=A{i};
            A{i}=A{j};
            A{j}=temp;
            %Record the swap in the index list.
            temp=idxList(i);
            idxList(i)=idxList(j);
            idxList(j)=temp;
        end
    end

    %Swap A{i+1} and A{r}
    temp=A{i+1};
    A{i+1}=A{r};
    A{r}=temp;
    
    %Record the swap in the index list.
    temp=idxList(i+1);
    idxList(i+1)=idxList(r);
    idxList(r)=temp;
    
    %The new pivot point to return.
    q=i+1;
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
