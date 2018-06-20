classdef queue < handle
%%QUEUE   A queue class. it manages a linked list implemented using the
%         dlNode class. Items from the top and bottom of the queue can be
%         directly accessed. 
%
%The concept of an queue is described in many introductory computer science
%texts, such as Chapter 5.4 of [1]. The implementation is generic and uses
%a linked list. In some instances, the more efficient but somewhat more
%limited priority_queue class (discussed therein) is preferable to use.
%
%Note that in many applications, it is better to just preallocate an array
%to use as a queue and to use an index to identify the current element in
%the array that is at the top of the queue. A queue like this one will
%cause memory allocation and deallocation events every time something is
%added or deleted, which might be slow.
%
%REFERENCES:
%[1] W. Ford, and W. Topp, Data Structures with C++. Upper Saddle River,
%    NJ: Prentice Hall, 1996.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        numInQueue
        top
        bottom
    end
    methods
       
      function newQueue = queue(initData)
      %QUEUE Create a new instance of the queue class.
      %
      %INPUTS: initData  Optional data that can form the first element in
      %                  the queue.
      %        
      %OUTPUTS: newQueue A new instance of the queue class that contains
      %                  initData, if supplied. Otherwise, newQueue is am
      %                  empty queue.
 
        if nargin > 0
            newQueue.top= dlnode(initData);
            newQueue.bottom=newQueue.top;
            newQueue.numInQueue=1;
        else
            newQueue.top=[];
            newQueue.bottom=[];
            newQueue.numInQueue=0;
        end
      end
      
      function boolVal=isEmpty(curQueue)
      %%ISEMPTY Return a boolean value indicating whether there is nothing
      %         in the queue.
      
          boolVal=(curQueue.numInQueue==0);
      end
      
      function newQueue = duplicate(curQueue)
      %DUPLICATE Duplicate the queue.
      %
      %        
      %OUTPUTS: newQueue A new instance of the queue class that contains
      %                  the same data as the calling queue.
      %
      %This function does not call a duplicate function on the data in the
      %queue, so if the queue is full of pointers (handle class
      %instances), the data in the new queue is the same set of pointers.
          
          if(curQueue.numInQueue==0)%If this queue is empty
              newQueue=queue();
          else%If there is data in the queue, then copy all of it into the
              %new queue.
              newQueue=queue(curQueue.bottom.Data);
              curItem=curQueue.bottom.Next;
              while(~isempty(curItem))
                 newQueue.addToTop(curItem.Data);
                 curItem=curItem.Next;
              end
          end
      end
      
      function theArray= cellArray(curQueue)
      %%CELLARRAY  Copy all of the data in the queue into a cell array.
      %            Pointers (handle class instances) in the queue are just
      %            copied; the underlying data is not duplicated.
          
          %This function condenses all of the data from the linked list into a cell array. 
          theArray=cell(curQueue.numInQueue,1);
          curItem=curQueue.bottom;
          curStep=1;
          while ~isempty(curItem)
            theArray{curStep}=curItem.Data;
            curStep=curStep+1;
            curItem=curItem.Next;
          end
      end
       
      function addToTop(curQueue,newData)
      %%ADDTOTOP Add the given data to the top of the queue.
      
        newNode=dlnode(newData);
        if(curQueue.numInQueue>0)
            newNode.insertAfter(curQueue.top);
            curQueue.top=newNode;
        else
            curQueue.top=newNode;
            curQueue.bottom=newNode;
        end
        curQueue.numInQueue=curQueue.numInQueue+1;
      end
      
      function addToBottom(curQueue,newData)
      %%ADDTOTOP Add the given data to the bottom of the queue.
      
        newNode=dlnode(newData);
        if(curQueue.numInQueue>0)
            newNode.insertBefore(curQueue.bottom);
            curQueue.bottom=newNode;
        else
            curQueue.top=newNode;
            curQueue.bottom=newNode;
        end
        curQueue.numInQueue=curQueue.numInQueue+1;
      end
      
      function data=getTopNode(curQueue)
      %%GETTOPNODE  Return the top dlnode instance in the queue (not just
      %             the data in the top node) without removing it from the
      %             queue.
      
          if(curQueue.numInQueue>0)
              data=curQueue.top;
          else
              data=[];
          end
      end
      
      function data=getBottomNode(curQueue)
      %%GETBOTTOMNODE  Return the bottom dlnode instance in the queue (not
      %                just the data in the bottom node) without removing
      %                it from the queue.
      
          if(curQueue.numInQueue>0)
              data=curQueue.bottom;
          else
              data=[];
          end
      end
      
      function data=dequeueBottom(curQueue)
      %%DEQUEUEBOTTOM Remove the bottom item from the queue and return the
      %               item.
      
        if(curQueue.numInQueue==0)%If there is nothing in the queue.
            data=[];
        else
            oldBottom=curQueue.bottom;
            curQueue.bottom=oldBottom.Next;
            data=oldBottom;
            oldBottom.disconnect;
            curQueue.numInQueue=curQueue.numInQueue-1;
            if(curQueue.numInQueue==0)
                curQueue.top=[];
            end
        end
      end
      
      function data=dequeueTop(curQueue)
      %%DEQUEUETOP Remove the top item from the queue and return the item.
          
        if(curQueue.numInQueue==0)
            data=[];
        else
            oldTop=curQueue.top;
            curQueue.top=oldTop.Prev;
            data=oldTop;
            oldTop.disconnect;
            curQueue.numInQueue=curQueue.numInQueue-1;
            if(curQueue.numInQueue==0)
                curQueue.bottom=[];
            end
        end
      end
      
      function deleteQueueItem(curQueue,item)
      %%DELETEQUEUEITEM Given a pointer to a dlnode in the queue, remove
      %                 that node from the queue and delete it. This does
      %                 not delete the data within the node.
      
          curQueue.numInQueue=curQueue.numInQueue-1;
          if(curQueue.bottom==item)
             curQueue.bottom=item.Next;
          end

          if(curQueue.top==item)
              curQueue.top=item.Prev;
          end
          
          item.delete;
      end
      
      function disp(curQueue)
      %DISP Display information about the queue.
          
        if(curQueue.numInQueue==0)
            disp('The queue is empty.');
        else
            disp('A queue containing:')
            curQueue.numInQueue
            disp('items with top node:')
            curQueue.top.disp();
        end
      end
      
      function delete(curQueue)
      %%DELETE  Delete the queue. This gets rid of all items in the queue,
      %         but does not delete the data within the items if they are
      %         from the handle subclass.
      
          while(curQueue.numInQueue~=0)
              curQueue.dequeueBottom;
          end
      end
      
      function deleteWithData(curQueue)
      %%DELETEWITHDATA  Delete the queue and delete all of the items within
      %                 the queue assuming that they are descended from a
      %                 handle subclass.

          while(curQueue.numInQueue~=0)
              Data=curQueue.dequeueBottom;
              Data.delete;
          end
          curQueue.delete;
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
