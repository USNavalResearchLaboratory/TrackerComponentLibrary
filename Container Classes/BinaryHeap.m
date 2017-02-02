classdef BinaryHeap < handle
%%BINARYHEAP   An implementation of a priority queue as a binary heap. A
%              priority queue is a data structure that allows one to easily
%              get/ remove the maximum or minimum item in the queue.
%
%DEPENDENCIES: KeyVal.m
%
%The methods of the class are generally based on the implementation
%described in Chapter 6.4 of [1].
%
%The class makes use of the KeyValue class for storing keys associated with
%values. The entire heap is stored in an array, which can be preallocated
%to the maximum size and thus be efficient.
%
%REFERENCES:
%[1] M.A.Weiss, Data Structures and Algorithm Analysis in C++, 2nd ed.
%    Reading, MA: Addison-Wesley, 1999.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        numInHeap
        heapArray
        isMaxHeap
    end
    
    methods
        function newHeap=BinaryHeap(initialMaxSize,isMaxHeap)
        %BINARYHEAP Allocate space for a binary heap and specify whether it
        %           is a max heap.
        %
        %INPUT: initialMaxSize  The amount of room preallocated for the
        %                       binary heap. If omitted, then a size-1 heap
        %                       is preallocated.
        %       isMaxHeap       A boolean value indicating whether the top
        %                       of the heap is a maximum or a minimum
        %                       value. The default if not specified is
        %                       maximum.
        %
        %OUTPUTS: newHeap       The newly created heap.
        
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
        end
        
        function val=heapSize(curHeap)
        %HEAPSIZE  Return the number of elements in the heap.
           val=curHeap.numInHeap; 
        end
        
        function val=isEmpty(curHeap)
        %ISEMPTY  Return true if there are no nodes in the heap.
           val=curHeap.numInHeap==0; 
        end
        
        function buildHeapFromKeysData(curHeap,keyArray,dataArray)
        %BUILDHEAPFROMDATA  Build the heap from given data. If the heap was
        %                   already initialized, it will be destroyed.
        %
        %INPUTS:    curHeap       The implicitly passed object.
        %           keyArray      A kX1 array of key values.
        %           dataArray     A linear cell array or a vector
        %                         containing the data associated with each
        %                         key value.
            
            numKeys=length(keyArray);
            curHeap.numInHeap=numKeys;
            curHeap.heapArray(numKeys,1)=KeyVal;
            
            if(isa(dataArray,'cell'))
                for curKey=1:numKeys
                    curHeap.heapArray(curKey)=KeyVal(keyArray(curKey),dataArray{curKey});
                end
            else
                for curKey=1:numKeys
                    curHeap.heapArray(curKey)=KeyVal(keyArray(curKey),dataArray(curKey));
                end
            end
            
            idx=floor(curHeap.numInHeap/2);

            while(idx>0)
                curHeap.percolateDown(idx);
                idx=idx-1; 
            end
        end
        
        function insert(curHeap,key,value)
        %INSERT Insert an entry with a specific key and value into the heap.  
          
            curHeap.numInHeap=curHeap.numInHeap+1;
            hole=curHeap.numInHeap;
            
            if(curHeap.isMaxHeap)
                while(hole>1&&key>curHeap.heapArray(floor(hole/2)))
                    curHeap.heapArray(hole)=curHeap.heapArray(floor(hole/2));
                    hole=floor(hole/2);
                end
            else
                while(hole>1&&key<curHeap.heapArray(floor(hole/2)))
                    curHeap.heapArray(hole)=curHeap.heapArray(floor(hole/2));
                    hole=floor(hole/2);
                end
            end
            curHeap.heapArray(hole)=KeyVal(key,value);
        end
        
        function val=getTop(curHeap)
        %GETTOP Return the key-value pair at the top of the heap as a
        %       KeyVal object. If the heap is a max heap, then it is the
        %       element with the largest key. Otherwise, it is the element
        %       with the smallest key. If the heap is empty, then an empty
        %       matrix is returned.
        
            if(curHeap.numInHeap>0)
                val=curHeap.heapArray(1);
            else
                val=[];
            end
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
                    else
                        break;
                    end

                    hole=child; 
                end
            end
            
            curHeap.heapArray(hole)=temp;
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
