classdef dlnode < handle
%%DLNODE  A class to represent a doubly-linked list node. Multiple dlnode
%         objects may be linked together to create linked lists. Each node
%         contains a piece of data and provides access to the next and
%         previous nodes. The basis for this class is taken from the Matlab
%         documentation. Added to the code from the documentation are
%         comparison routines.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

   properties
      Data
      Next
      Prev
   end
    
   methods
        function node = dlnode(Data)
        % DLNODE  Constructs a dlnode object.
            if nargin > 0
                node.Data = Data;
            end
            node.Next=[];
            node.Prev=[];
        end
      
%The following funtions are so that multiple nodes can be compared to each
%other according to their data and directly to a data object.
        function val=lt(obj1,obj2)
        %%LT Less than comparison  
          
            if(isa(obj1,'dlnode')==true&&isa(obj2,'dlnode')==true)
                val=obj1.Data<obj2.Data;
            elseif(isa(obj1,'dlnode')==true)
                val=obj1.Data<obj2;
            else
                val=obj1<obj2.Data;
            end
        end
      
        function val=gt(obj1,obj2)
        %%GT Greater than comparison
        
            if(isa(obj1,'dlnode')==true&&isa(obj2,'dlnode')==true)
                val=obj1.Data>obj2.Data;
            elseif(isa(obj1,'dlnode')==true)
                val=obj1.Data>obj2;
            else
                val=obj1>obj2.Data;
            end
        end
      
        function val=le(obj1,obj2)
        %%LE Less than or equal to comparison
        
            if(isa(obj1,'dlnode')==true&&isa(obj2,'dlnode')==true)
                val=obj1.Data<=obj2.Data;
            elseif(isa(obj1,'dlnode')==true)
                val=obj1.Data<=obj2;
            else
                val=obj1<=obj2.Data;
            end
        end
      
        function val=ge(obj1,obj2)
        %%GE Greater than or equal to comparison
        
            if(isa(obj1,'dlnode')==true&&isa(obj2,'dlnode')==true)
                val=obj1.Data>=obj2.Data;
            elseif(isa(obj1,'dlnode')==true)
                val=obj1.Data>=obj2;
            else
                val=obj1>=obj2.Data;
            end
        end
      
        function val=ne(obj1,obj2)
        %%NE Not equal to comparison
        
            if(isa(obj1,'dlnode')==true&&isa(obj2,'dlnode')==true)
                val=obj1.Data~=obj2.Data;
            elseif(isa(obj1,'dlnode')==true)
                val=obj1.Data~=obj2;
            else
                val=obj1~=obj2.Data;
            end
        end
      
        function val=eq(obj1,obj2)
        %%EQ Equal to comparison
        
            if(isa(obj1,'dlnode')==true&&isa(obj2,'dlnode')==true)
                val=obj1.Data==obj2.Data;
            elseif(isa(obj1,'dlnode')==true)
                val=obj1.Data==obj2;
            else
                val=obj1==obj2.Data;
            end
        end
      
      function insertAfter(newNode, nodeBefore)
      % INSERTAFTER  Inserts newNode after nodeBefore.
      
         disconnect(newNode);
         newNode.Next = nodeBefore.Next;
         newNode.Prev = nodeBefore;
         if ~isempty(nodeBefore.Next)
            nodeBefore.Next.Prev = newNode;
         end
         nodeBefore.Next = newNode;
      end
      
      function insertBefore(newNode, nodeAfter)
      % INSERTBEFORE  Inserts newNode before nodeAfter.
      
         disconnect(newNode);
         newNode.Next = nodeAfter;
         newNode.Prev = nodeAfter.Prev;
         if ~isempty(nodeAfter.Prev)
             nodeAfter.Prev.Next = newNode;
         end
         nodeAfter.Prev = newNode;
      end 

      function disconnect(node)
      %DISCONNECT  Removes a node from a linked list. The node can be
      %            reconnected or moved to a different list. This does not
      %            delete the data in the node.
      
         if(~isscalar(node))
            error('Nodes must be scalar')
         end
  
         prevNode = node.Prev;
         nextNode = node.Next;
         if(~isempty(prevNode))
             prevNode.Next = nextNode;
         end
         if(~isempty(nextNode))
             nextNode.Prev = prevNode;
         end
         node.Next = [];
         node.Prev = [];
      end
      
      function delete(node)
      % DELETE  Deletes a dlnode from a linked list.
      
         disconnect(node);
      end  
      
      function disp(node)
      % DISP  Display a link node.
      
         if (isscalar(node))
            disp('Doubly-linked list node with data:')
            disp(node.Data)
         else % If node is an object array, display dims only
            dims = size(node);
            ndims = length(dims);
            for k = ndims-1:-1:1
               dimcell{k} = [num2str(dims(k)) 'x'];
            end
            dimstr = [dimcell{:} num2str(dims(ndims))];
            disp([dimstr ' array of doubly-linked list nodes']);
         end
      end 
   end % methods
end % classdef

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
