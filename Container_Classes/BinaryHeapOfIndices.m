classdef BinaryHeapOfIndices < handle
%%BINARYHEAPOFINDICES A binary heap stores values associated with keys. At
%    any time, one can easily access either the largest or the smallest
%    key, depending on the type of binary heap. This class is specialized
%    for the case where the values are indices starting from 1 that can be
%    inserted in any order. Whereas, in a regular binary heap, there is no
%    good way to change the key of an element in the heap given its value,
%    in this class, the array idxList is maintained so that is possible.
%    idxList is the length of all of the possible values and idxList(i)
%    holds the position of the value index i in the heap. If nothing with
%    index value i is in the heap, then idxList(i) is zero. This type of
%    binary heap is useful for implementing Dijkstra's algorithm, where one
%    has to increment things by value, and one knows that all values from 1
%    to some maximum size will be present (in a fully connected graph).
%
%DEPENDENCIES: KeyVal.m
%
%The methods of the class are generally based on the implementation
%described in Chapter 6.4 of [1].
%
%The class makes use of the KeyValue class for storing keys associated with
%values. The entire heap is stored in an array, which can be preallocated
%to the maximum size and thus be efficient. The class has been modified to
%keep track of idxList. The heap can be built by repeatedly calling the
%insert method. The method changeIndexedKey is not in traditional heaps;
%it lets one can change the key associated with a particular index value.
%
%REFERENCES:
%[1] M.A.Weiss, Data Structures and Algorithm Analysis in C++, 2nd ed.
%    Reading, MA: Addison-Wesley, 1999.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        numInHeap
        heapArray
        isMaxHeap
        idxList
    end
    
    methods
        function newHeap=BinaryHeapOfIndices(initialMaxSize,isMaxHeap)
        %BINARYHEAPOFINDICES Allocate space for a binary heap and specify
        %                    whether it is a max heap.
        %
        %INPUT: initialMaxSize The amount of room preallocated for the
        %                      binary heap. If omitted, then a size-1 heap
        %                      is preallocated.
        %            isMaxHeap A boolean value indicating whether the top
        %                      of the heap is a maximum or a minimum
        %                      value. The default if not specified is
        %                      maximum.
        %
        %OUTPUTS: newHeapThe newly created heap.
        
            if(nargin<2)
                isMaxHeap=true;
            end
            if(nargin<1)
                initialMaxSize=1;
            end

            %This line is here, because the next line will fail if Matlab
            %has not identified heapArray as being an array of instances of
            %the KeyValue class.
            newHeap.heapArray=KeyVal();
            newHeap.heapArray(initialMaxSize,1)=KeyVal();
            newHeap.numInHeap=0;
            newHeap.isMaxHeap=isMaxHeap;
            newHeap.idxList=[];
        end
        
        function val=heapSize(curHeap)
        %HEAPSIZE  Return the number of elements in the heap.
           val=curHeap.numInHeap; 
        end
        
        function val=isEmpty(curHeap)
        %ISEMPTY  Return true if there are no nodes in the heap.
           val=curHeap.numInHeap==0; 
        end
        
        function boolVal=indexIsInHeap(curHeap,index)
        %INDEXISINHEAP Returns true if a key with a value having the given
        %              index is already in the heap.
            if(index>length(curHeap.idxList)||curHeap.idxList(index)==0)
                boolVal=false;
            else
                boolVal=true;
            end
        end
        
        function insert(curHeap,key,index2Add)
        %INSERT Insert an entry with a specific key and index value into 
        %       the heap.  
            
            %Insert the key-value pair into the heap.
            curHeap.numInHeap=curHeap.numInHeap+1;
            hole=curHeap.numInHeap;
            curHeap.heapArray(hole)=KeyVal(key,index2Add);
            curHeap.idxList(index2Add)=hole;
            
            %Adjust the heap to put the key-value pair into the correct
            %position. 
            curHeap.percolateUp(hole);
        end
        
        function val=getTop(curHeap)
        %GETTOP Return the key-value pair at the top of the heap as a
        %       KeyVal object. If the heap is a max heap, then it is the
        %       element with the largest key. Otherwise, it is the element
        %       with the smallest key. If the heap is empty, then an empty
        %       matrix is returned.
        %
        %The value returned is a shallow copy of the KeyVal object in the
        %heap. Changes to it do not affect the values in the heap, unless
        %those values are handle class objects. Do not try modifying the
        %key of an object in the heap by modifying the returned KeyVal
        %object. Use the function changeIndexedKey for changing keys.
        
            if(curHeap.numInHeap>0)
                val=copy(curHeap.heapArray(1));
            else
                val=[];
            end
        end
        
        function didSucceed=changeIndexedKey(curHeap,newKey,index)
        %%CHANGEINDEXEDKEY Change the key of an item with a particular
        %                  index value in the heap and then adjust the
        %                  ordering of the heap. If the item is not in the
        %                  heap, then didSucceed will be returned false,
        %                  otherwise it will be true. 
            
            %If the indexed item is not in the heap.
            if(~curHeap.indexIsInHeap(index))
                didSucceed=false;
                return;
            end
            
            %Otherwise, find the item in the heap.
            idx2Change=curHeap.idxList(index);
            
            oldKey=curHeap.heapArray(idx2Change).key;
            curHeap.heapArray(idx2Change).key=newKey;
            
            %We now have to rebalance the heap. If the element may have
            %become too large for its position in a min heap or too small
            %for its position in a max heap, then it has to be percolated
            %down. otherwise, it might have to be percolated up.
            if((newKey>oldKey&&curHeap.isMinHeap)||(newKey<oldKey&&curHeap.isMaxHeap))
                curHeap.percolateDown(idx2Change);
            else%Otherwise, it has to be percolated up.
                curHeap.percolateUp(idx2Change);
            end
            
            didSucceed=true;
        end
        
        function val=deleteTop(curHeap)
        %DELETETOP Remove the top element of heap heap (the one that getTop
        %          would return) and return the key and value of the
        %          element as a KeyVal object.
            
            if(curHeap.numInHeap==0)
                val=[];
                return;
            end
            val=curHeap.heapArray(1);
            
            curHeap.heapArray(1)=curHeap.heapArray(curHeap.numInHeap);
            curHeap.numInHeap=curHeap.numInHeap-1;
            curHeap.percolateDown(1);
            
            %Mark the index as not being in the heap.
            removedIdx=val.value;
            curHeap.idxList(removedIdx)=0;
        end
    end
    
    methods(Access=private)        
        function percolateDown(curHeap,hole)
            temp=curHeap.heapArray(hole);
            
            if(curHeap.isMaxHeap)
                while(2*hole<=curHeap.numInHeap)
                    child=2*hole;

                    if(child~=curHeap.numInHeap&&curHeap.heapArray(child+1)>curHeap.heapArray(child))
                       child = child+1; 
                    end

                    if(curHeap.heapArray(child)>temp)
                        curHeap.heapArray(hole)=curHeap.heapArray(child);
                        
                        %Update the indexed position of the item in the
                        %heap.
                        movedIdx=curHeap.heapArray(hole).value;
                        curHeap.idxList(movedIdx)=hole;
                    else
                        break;
                    end

                    hole=child; 
                end
            else
                while(2*hole<=curHeap.numInHeap)
                    child=2*hole;

                    if(child~=curHeap.numInHeap&&curHeap.heapArray(child+1)<curHeap.heapArray(child))
                       child = child+1; 
                    end

                    if(curHeap.heapArray(child)<temp)
                        curHeap.heapArray(hole)=curHeap.heapArray(child);
                        %Update the indexed position of the item in the
                        %heap.
                        movedIdx=curHeap.heapArray(hole).value;
                        curHeap.idxList(movedIdx)=hole;
                    else
                        break;
                    end

                    hole=child; 
                end
            end
            
            curHeap.heapArray(hole)=temp;
            
            %Update the indexed position of the item in the heap.
            movedIdx=curHeap.heapArray(hole).value;
            curHeap.idxList(movedIdx)=hole;
        end
    
        function percolateUp(curHeap,hole)
            temp=curHeap.heapArray(hole);
            if(curHeap.isMaxHeap)
                floorIdx=floor(hole/2);
                while(hole>1&&temp.key>curHeap.heapArray(floorIdx))
                    curHeap.heapArray(hole)=curHeap.heapArray(floorIdx);
                    %Update the indexed position of the item in the heap.
                    movedIdx=curHeap.heapArray(hole).value;
                    curHeap.idxList(movedIdx)=hole;
                    
                    hole=floorIdx;
                    floorIdx=floor(hole/2);
                end
            else
                floorIdx=floor(hole/2);
                while(hole>1&&temp.key<curHeap.heapArray(floorIdx))
                    curHeap.heapArray(hole)=curHeap.heapArray(floorIdx);
                    %Update the indexed position of the item in the heap.
                    movedIdx=curHeap.heapArray(hole).value;
                    curHeap.idxList(movedIdx)=hole;
                    
                    hole=floorIdx;
                    floorIdx=floor(hole/2);
                end
            end
            curHeap.heapArray(hole)=temp;
            
            %Update the indexed position of the item in the heap.
            movedIdx=temp.value;
            curHeap.idxList(movedIdx)=hole;
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
