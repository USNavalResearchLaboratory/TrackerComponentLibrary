classdef CircularlyLinkedList < handle
%%CIRCULARLYLINKEDLIST A class to hold a doubly circularly linked list. The
%                    next node after the last node is thus the first node.
%                    Each node in the list can hold some type of data. The
%                    class keeps a pointer to the current node in
%                    the list, which can be moved forward or backward.
%                    Methods can be used to insert data before or after the
%                    current node and to delete the current node, which
%                    means the other nodes will be connected to each other.
%
%The concepts behind linked lists and circularly linked lists are covered
%in many introductory computer science textbooks.
%
%When the linked list is deleted, delete is not explicitely called for the
%data within the nodes, though Matlab's garbage collection should
%eventually call the destructors if nothing else points to them.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        numNodesInList
        currentNode
    end
    methods
        function newList=CircularlyLinkedList(initialData)
        %%CIRCULARLINKEDLIST The constructor method of the
        %                    CircularlyLinkedList class.
        %
        %INPUTS: initialData Data for an initial node to be put into the
        %                    list. If this is omitted, an empty list is
        %                    created.
        %
        %OUTPUS:newList  A new CircularlyLinkedList object.
        
            if(nargin<1)
                newList.numNodesInList=0;
                newList.currentNode=[];
            else
                newList.numNodesInList=1;
                newList.currentNode=dlnode(initialData);
                newList.currentNode.Next=newList.currentNode;
                newList.currentNode.Prev=newList.currentNode;
            end
        end
        
        function val=listLength(curList)
        %LISTLENGTH  Return the number of elements in the list.
           val=curList.numNodesInList; 
        end
        
        function val=isEmpty(curList)
        %ISEMPTY  Return true if there are no nodes in the list.
            val=curList.numNodesInList==0;
        end
        
        function insertAfter(curList,theData)
        %%INSERTAFTER Insert the data in a new node after the current node
        %             in the list. The current node is not changed.
            
            if(curList.numNodesInList==0)
                curList.currentNode=dlnode(theData);
                curList.currentNode.Next=curList.currentNode;
                curList.currentNode.Prev=curList.currentNode;
                curList.numNodesInList=1;
            else
                nextNode=curList.currentNode.Next;
                newNode=dlnode(theData);
                newNode.Prev=curList.currentNode;
                newNode.Next=nextNode;
                curList.currentNode.Next=newNode;
                nextNode.Prev=newNode;

                curList.numNodesInList=curList.numNodesInList+1;
            end
        end
        
        function insertBefore(curList,theData)
        %%INSERTAFTER Insert the data in a new node before the current node
        %             in the list. The current node is not changed.
            
            if(curList.numNodesInList==0)
                curList.currentNode=dlnode(theData);
                curList.currentNode.Next=curList.currentNode;
                curList.currentNode.Prev=curList.currentNode;
                curList.numNodesInList=1;
            else
                prevNode=curList.currentNode.Prev;
                newNode=dlnode(theData);
                curList.currentNode.Prev=newNode;
                newNode.Next=curList.currentNode;
                newNode.Prev=prevNode;
                prevNode.Next=newNode;
                
                curList.numNodesInList=curList.numNodesInList+1;
            end
        end
        
        function theData=getCurData(curList)
        %%GETCURDATA Get the data in the current node. If the list is
        %            empty, an empty matrix is returned.
            if(curList.numNodesInList==0)
                theData=[];
            else
                theData=curList.currentNode.Data;
            end
        end
        
        function goToNextNode(curList)
        %%GOTONEXTNODE Change the current node to the next node in the
        %              circularly linked list. If the list is empty, there
        %              will be an error. 
            curList.currentNode=curList.currentNode.Next;
        end
        
        function goToPrevNode(curList)
        %%GOTOPREVNODE Change the current node to the previous node in the
        %              circularly linked list. If the list is empty, there
        %              will be an error. 
            
            curList.currentNode=curList.currentNode.Prev;
        end
        
        function deleteCurrentNode(curList)
        %%DELETECURRENTNODE Delete the current node in the list. The new
        %                   current node after deletion will be what was
        %                   the next node. If this list is empty, then
        %                   there will be an error.
        
            curList.numNodesInList=curList.numNodesInList-1;
            if(curList.numNodesInList>0)
                curList.currentNode.Prev.Next=curList.currentNode.Next;
                curList.currentNode.Next.Prev=curList.currentNode.Prev;
                temp=curList.currentNode.Next;
                curList.currentNode.delete();
                curList.currentNode=temp;
            else
                curList.currentNode=[];
            end
        end
        
        function delete(curList)
            %The destructor
            while(curList.numNodesInList>0)
                curList.deleteCurrentNode();
            end
        end
        
        function disp(curList)
        %DISP Display information about the list.
          
        if(curList.numNodesInList==0)
            disp('The circularly linked list is empty.');
        else
            disp('A circularly linked list containing:')
            curList.numNodesInList
            disp('current node data:')
            curList.currentNode.Data
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
