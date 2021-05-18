classdef AVLTree < handle
%%AVLTREE A class for an Adelson-Velskii, Landis (AVL) binary search tree.
%         This tree remains balanced in height after each insertion or
%         deletion. The tree stores key-value pairs that one can rapidly
%         search. Duplicate keys are not permitted in the tree.
%
%The basic implementation is somewhat based on the description in Chapter
%10 of [1]. However, AVL trees are described in most introductory books on
%data structures. Unlike the implementation described in Preiss, this class
%is not implemented as a subclass of more general binary tree data
%structures.Also, no needless empty classes are generated to mark the
%leaves of the tree. This gets rid of the need for special Left and Right
%accessor methods.
%
%The AVL tree stores KeyVal objects and allows one to search for them by
%key. In a practical (i.e. C/C++) implementation where an upper bound on
%the number of KeyVal pairs in the tree can be set, it would make sense to
%pair the AVLTree class with a custom memory allocation routine for all of
%the entries in the tree. That is, one allocated all of the possible nodes
%in the tree as a single block of memory in advance and puts all of the
%pointers to the nodes in a binary heap (allocated as a single array). Each
%time a node is needed for the tree, it is popped from the heap. Each time
%a node is freed, it is pushed back on the heap. This should, in theory, be
%faster than calling malloc and free many times, as they normally maintain
%more complicated memory management structures.
%
%A C++ implementation of this class for Matlab is difficult to make, since
%Matlab allows for more than just pointers or scalars as the values in the
%KeyVal class. Due to the overhead of classes in Matlab, it is not
%guaranteed that using an AVL tree is going to be faster than if one put
%the key-value pairs in an array and resorted the array every single time
%something is added. However, when dealing with a large amount of data,
%that is not necessarily the smartest approach when one wants a dynamic
%searchable data structure in C/C++.
%
%REFERENCES:
%[1] B. R. Preiss, Data Structures and Algorithms with Object-Oriented
%    Design Patterns in C++. New York, NY: John Wiley & Sons, Inc., 1999.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
    %Duplicate keys are not permitted
    properties
        height
        keyVal
        left
        right
    end
    
    methods
        function theTree=AVLTree(keyValInit)
        %%AVLTREE The constructor method of the AVLTree class.
        %
        %INPUTS: keyValInit An instance of the keyVal class that holds the
        %                   key and data for the new tree/ new node. If
        %                   this parameter is omitted, then an empty tree
        %                   is created.
        %
        %OUTPUTS: theTree  A new AVLTree object.
        %
        %The initialization routine is based on page 319 of the Preiss
        %book.
            
            %The constructor
            if(nargin<1)
                theTree.keyVal=[];
                theTree.height=-1;
            else
                theTree.keyVal=keyValInit;
                theTree.height=0;
            end
            
            theTree.left=[];
            theTree.right=[];
        end
        
        function theTree=clear(theTree)
        %%CLEAR Clear the tree.
        %
        %
        %OUTPUTS: theTree  An empty AVLTree object.
        %
        %This function clears all the KeyVals stored in the tree.
            
            theTree.keyVal=[];
            theTree.height=-1;
            theTree.left=[];
            theTree.right=[];
        end

        function val=isempty(theTree)
        %%ISEMPTY This returns true if there are no KeyVal pairs in this
        %         tree.
            
            val=(theTree.height==-1);
        end
        
        %Accessor methods
        function val=KeyVal(theTree)
        %%KEYVAL Get the KeyVal object with which the instance was
        %        initialized. This is empty if there is nothing in the
        %        tree. If nothing is in the tree, then this is an empty
        %        matrix.
        
            val=theTree.keyVal;
        end
        
        function insert(theTree,keyVal2Insert)
        %%INSERT Insert a new keyVal pair into the tree. An error is
        %        raised if a node with the same key already exists in the
        %        tree.
        %
        %INPUTS: theTree The implicitly passed AVLTree object.
        %     key2Insert An instance of the KeyVal class containing the
        %                key and its associated value to insert.
        %
        %OUTPUTS: none
        %
        %The KeyVal object is directly inserted into the tree without
        %duplicating it. Subsequent changes to the passed KeyVal object
        %will affect what is in the tree. Changing the value of the object
        %will invalidate the tree structure.
        %
        %This function is based on page 312 of the Preiss book.
            
            if(isempty(theTree.keyVal))
                %If this node is empty, then just attach the key.
                theTree.height=0;
                theTree.keyVal=keyVal2Insert;
                return;
            else
                if(keyVal2Insert<theTree.keyVal)
                    if(isempty(theTree.left))
                        %Attach it to the empty left node
                        theTree.left=AVLTree(keyVal2Insert);
                    else
                        theTree.left.insert(keyVal2Insert);
                    end
                    
                elseif(keyVal2Insert>theTree.keyVal)
                    if(isempty(theTree.right))
                        %Attach the key to the empty right node
                        theTree.right=AVLTree(keyVal2Insert);
                    else
                        theTree.right.insert(keyVal2Insert);
                    end
                else
                    error('Attempted to insert a duplicate key!')
                end
            end
            
            theTree.balance();
        end
        
        function didSucceed=remove(theTree,theKey)
        %%REMOVE Remove a node having the given key from the tree.
        %
        %INPUTS: theTree The implicitly passed AVLTree object.
        %         theKey The key of the object to remove from the tree.
        %
        %OUTPUTS: didSucceed A boolean value indicating whether the key was
        %                    removed from the tree. This will only be false
        %                    if the key was not in the tree.
        %
        %This is based on the description on page 314 of Preiss, except it
        %allows for a failure to find the node. Also, unlike in Preiss,
        %nodes are completely removed from the tree rather than left empty,
        %except for a node at the root of the tree, which will just be made
        %empty if this function is called.
            
            %If the tree is empty, then there is nothing to get rid of.
            if(isempty(theTree.keyVal))
                didSucceed=false;
                return;
            end
            
            if(theKey==theTree.keyVal.key)%If this node should be removed.
                if(~isempty(theTree.left))
                    theTree.keyVal=theTree.left.findMax();
                    theTree.left.remove(theTree.keyVal);
                    
                    %If the left node should be deleted.
                    if(isempty(theTree.left.keyVal))
                        delete(theTree.left);
                        theTree.left=[];
                    end
                elseif(~isempty(theTree.right))
                    theTree.keyVal=theTree.right.findMin();
                    theTree.right.remove(theTree.keyVal);
                    
                    %If the right node should be deleted.
                    if(isempty(theTree.right.keyVal))
                        delete(theTree.right);
                        theTree.right=[];
                    end
                else
                    theTree.keyVal=[];
                    theTree.height=-1;
                    didSucceed=true;
                    return;
                end
                didSucceed=true;
            elseif(theKey<theTree.keyVal.key&&~isempty(theTree.left))
                didSucceed=theTree.left.remove(theKey);
                
                %If the left node should be deleted.
                if(isempty(theTree.left.keyVal))
                    delete(theTree.left);
                    theTree.left=[];
                end
            elseif(~isempty(theTree.right))
                didSucceed=theTree.right.remove(theKey);
                
                %If the right node should be deleted.
                if(isempty(theTree.right.keyVal))
                    delete(theTree.right);
                    theTree.right=[];
                end
            else
                didSucceed=false;
                return;
            end
            
            %Rebalance the tree.
            theTree.balance();
        end
        
        function foundKeyVal=find(theTree,key2Find)
        %%FIND Find the KeyVal object corresponding to a key in the tree.
        %      If a key is not in the tree, then an empty matrix is
        %      returned.
        %
        %INPUTS: theTree The implicitly passed AVLTree object.
        %       key2Find The key to find in the tree.
        %
        %OUTPUTS: foundKeyVal The keyVal object having the provided key. If
        %                     the key is not in the tree, then an empty
        %                     matrix is returned.
        %
        %The returned KeyVal class is the actual value in the tree, which
        %is a handle class. The key itself should NOT be modified. Doing so
        %would invalidate the structure of the AVLTree. If one wishes to
        %modify a key, they should remove the old key/value from the tree
        %and then reinsert it under a new key.
        %
        %This is based on the description on page 311 of the Preiss book.
            
            %If there is nothing in the tree.
            if(isempty(theTree.keyVal))
                foundKeyVal=[];
                return;
            end
            
            %If the key matches this node, return the value of this
            %node.
            if(key2Find==theTree.keyVal.key)
                foundKeyVal=theTree.keyVal;
            elseif(key2Find<theTree.keyVal.key&&~isempty(theTree.left))
                foundKeyVal=theTree.left.find(key2Find);
            elseif(~isempty(theTree.right))
                foundKeyVal=theTree.right.find(key2Find);
            else
                %There are no left or right children and this node is not
                %the key, so the key is not in the tree.
                foundKeyVal=[];
            end
        end
        
        function keyNotFound=replaceValue4Key(theTree,key2Find,newValue)
        %%REPLACEVALUE4KEY Find the KeyVal object corresponding to a key in
        %       the tree and replace its value with the new given value
        %       (The key remains unchanged).
        %
        %INPUTS: theTree The implicitly passed AVLTree object.
        %       key2Find The key to find in the tree.
        %       newValue The new value to associate with the specified key.
        %
        %OUTPUTS: keyNotFound If the key is not found (and thus the tree
        %                remains unchanged, then this is true. Otherwise,
        %                this is false.
        %
        %This function just calls "find" and then changes the value part of
        %the KeyVal pair.
            
            foundKeyVal=theTree.find(key2Find);
            
            if(isempty(foundKeyVal))
                keyNotFound=true;
            else
                foundKeyVal.value=newValue;
                keyNotFound=false;
            end
        end
        
        function foundKeyVal=findMin(theTree)
        %%FINDMIN Find the minimum value in the tree. If the tree is empty,
        %         then an empty matrix will be returned.
        %
        %This is based on the description on page 311 of the Preiss book.
            
            %If there is nothing in the tree.
            if(isempty(theTree.keyVal))
                foundKeyVal=[];
                return;
            end
            
            if(isempty(theTree.left))
                foundKeyVal=theTree.keyVal;
            else
                foundKeyVal=theTree.left.findMin();
            end
        end
        
        function foundKeyVal=findMax(theTree)
        %FINDMAX Find the maximum value in the tree. If the tree is empty,
        %        then an empty matrix will be returned.
        %
        %This logically follows from how the foundMin function is
        %implemented.
            
            %If there is nothing in the tree.
            if(isempty(theTree.keyVal))
                foundKeyVal=[];
                return;
            end
            
            if(isempty(theTree.right))
                foundKeyVal=theTree.keyVal;
            else
                foundKeyVal=theTree.right.findMax();
            end
        end
        
        function inOrderTraversal(theTree,evalFuncHandle)
        %%INORDERTRAVERSAL A recursive function that visits all of
        %                  the nodes in the tree in order of increasing
        %                  keys. At each node, the function evalFuncHandle
        %                  is called with the keyVal pair. An in-order
        %                  traversal of an AVL tree is the same as a
        %                  depth-first traversal.
        %
        %INPUTS: theTree The implicitly passed AVLTree object. If the tree
        %                is empty, then evalFuncHandle will never be
        %                called.
        % evalFuncHandle A handle for a function that is called at every
        %                node with the keyVal object that forms the data of
        %                the node.
        %
        %OUTPUTS: none
        %
        %Note that in returning information via the evalFuncHandle, the
        %HandleWrapper class might be useful in allowing different levels
        %of recursion modify the same data without returning anything.
        %
        %This is based on the description of such a traversal from page 289
        %of the Preiss book.
            
            if(~isempty(theTree.left))
                theTree.left.inOrderTraversal(evalFuncHandle);
            end
            
            if(isempty(theTree.keyVal))
                return;
            end
            
            evalFuncHandle(theTree.keyVal);
            
            if(~isempty(theTree.right))
                theTree.right.inOrderTraversal(evalFuncHandle);
            end
        end
        
        function theKeys=getAllKeys(theTree)
        %%GETALLKEYS This function returns a cell array containing all of
        %            the keys below this node in the tree.
            
            %We do not know the total number of keys, but this is an upper
            %bound.
            maxNumKeys=2^(theTree.height+1)-1;
            
            keyList=[];
            keyList.numInList=0;
            keyList.theList=cell(maxNumKeys,1);
            %We use a handleWrapper so that the same data can be modified
            %in each of the recursive steps.
            keyListWrapper=HandleWrapper(keyList);
            evalFuncHandle=@(theKeyVal)addKey(theKeyVal,keyListWrapper);
            
            %Fill in the values.
            theTree.inOrderTraversal(evalFuncHandle);
            
            %Shrink to fit for return.
            theKeys=keyListWrapper.data.theList(1:keyListWrapper.data.numInList);
            
            function addKey(theKeyVal,theWrapper)
                theWrapper.data.numInList=theWrapper.data.numInList+1;
                theWrapper.data.theList{theWrapper.data.numInList}=theKeyVal.key;
            end
        end
        
        function val=count(theTree)
        %%COUNT Return the total number of KeyVal values in the tree.
            
            %IF the tree is empty
            if(theTree.height==-1)
                val=0;
                return
            end
            
            %We have to use a handle wrapper so that all of the recursions
            %in the inOrderTraversal function edit the same data.
            theWrapper=HandleWrapper(0);
            evalFuncHandle=@(x)incVal(x,theWrapper);
            
            %Count how many things are in the tree.
            theTree.inOrderTraversal(evalFuncHandle);
            val=theWrapper.data;
            
            function incVal(~,theWrapper)
                theWrapper.data=theWrapper.data+1;
            end
        end
        
        function delete(theTree)
        %DELETE The destructor function.
            
            if(~isempty(theTree.left))
                delete(theTree.left);
            end
            
            if(~isempty(theTree.right))
                delete(theTree.right);
            end
        end
    end
    
    
    methods(Access=private)
        function val=balanceFactor(theTree)
        %%BALANCEFACTOR Return how well balanced this node is. This is used
        %               to determine whether the node needs to be
        %               rebalanced.
        %
        %This is based on the description on page 319 of the Preiss book.
            
            %If there is nothing in the tree.
            if(isempty(theTree.keyVal))
                val=0;
                return;
            end
            
            if(isempty(theTree.left))
                leftHeight=-1;
            else
                leftHeight=theTree.left.height;
            end
            
            if(isempty(theTree.right))
                rightHeight=-1;
            else
                rightHeight=theTree.right.height;
            end
            val=leftHeight-rightHeight;
        end
        
        function adjustHeight(theTree)
        %%ADJUSTHEIGHT This function corrects the height of a node after
        %              rotations for rebalancing.
        %
        %This is based on the description on page 319 of the Preiss book.
            
            if(isempty(theTree.keyVal))
                theTree.height=-1;
            else
                if(isempty(theTree.left))
                    leftHeight=-1;
                else
                    leftHeight=theTree.left.height;
                end
                
                if(isempty(theTree.right))
                    rightHeight=-1;
                else
                    rightHeight=theTree.right.height;
                end
                
                theTree.height=1+max(leftHeight,rightHeight);
            end
        end
        
        function LLRotation(theTree)
        %%LLROTATION Perform an LL rotation of the node and its children.
        %            This assumes that the appropriate left child node is
        %            not empty.
        %
        %This is based on the description on page 324 of the Preiss book.
            
            temp=theTree.right;
            theTree.right=theTree.left;
            theTree.left=theTree.right.left;
            theTree.right.left=theTree.right.right;
            theTree.right.right=temp;
            
            tempKeyVal=theTree.keyVal;
            theTree.keyVal=theTree.right.keyVal;
            theTree.right.keyVal=tempKeyVal;
            
            theTree.right.adjustHeight();
            theTree.adjustHeight();
        end
        
        function RRRotation(theTree)
        %%RRROTATION Perform an RR rotation of the node and its children.
        %            This assumes that the appropriate right child node is
        %            not empty.
        %
        %This is just a mirror image of the LL rotation. That is, one
        %switches all left and right.
            
            temp=theTree.left;
            theTree.left=theTree.right;
            theTree.right=theTree.left.right;
            theTree.left.right=theTree.left.left;
            theTree.left.left=temp;
            
            tempKeyVal=theTree.keyVal;
            theTree.keyVal=theTree.left.keyVal;
            theTree.left.keyVal=tempKeyVal;
            
            theTree.left.adjustHeight();
            theTree.adjustHeight();
        end
        
        function LRRotation(theTree)
        %%LRROTATION Perform an LR rotation of the node and its children.
        %            This assumes that the appropriate child nodes are not
        %            empty.
        %
        %This is based on the description on page 325 of the Preiss book.
            
            theTree.left.RRRotation();
            theTree.LLRotation();
        end
        
        function RLRotation(theTree)
        %%RLROTATION Perform an RL rotation of the node and its children.
        %            This assumes that the appropriate child nodes are not
        %            empty.
        %
        %This is just a mirror image of the LR rotation. That is, one
        %switches all left and right.
            
            theTree.right.LLRotation();
            theTree.RRRotation();
        end
        
        function balance(theTree)
        %%BALANCE Balance this node in the tree using rotations, assuming
        %         that the lack of balance is no more than once node.
        %
        %This is based on the description on page 326 of the Preiss book.
            
            theTree.adjustHeight();
            
            balFactor=theTree.balanceFactor();
            
            if(abs(balFactor)>1)
                if(balFactor>0)
                    if(theTree.left.balanceFactor()>0)
                        theTree.LLRotation();
                    else
                        theTree.LRRotation();
                    end
                else
                    if(theTree.right.balanceFactor()<0)
                        theTree.RRRotation();
                    else
                        theTree.RLRotation();
                    end
                end
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
