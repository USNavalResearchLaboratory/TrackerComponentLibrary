classdef kdTree < handle
%%KDTREE A k-dimensional tree for performing orthogonal range queries and
%        nearest neighbor searches on a batch of data. If a C++ class
%        interface for the kd tree has been compiled, then its functions
%        will be called in place of the Matlab routines, since the C++
%        implementation can be significantly faster.
%
%The kd tree for range queries is described in [1]. The algorithm for a
%nearest neighbor search is described in [2]. The squared Euclidean
%distance is used as the metric for the nearest neighbor search.
%
%This implementation builds the tree from a batch of data all at once. The
%split point is chosen to be the median along the split dimension. Note
%that if the C++ implementation is used, the mex file is locked when a
%kdtree object is created and is not unlocked (and able to be recompiled)
%until all of the kdTree objects have been freed.
%
%Note that the datapoints are just stored in the data member of the class
%and are not reordered. Thus, the datapoints can remain associated with any
%additional information that is not passed to the class, but which is
%indexed the same way as the data points.
%
%The index of the root data point of the tree is stored at DATAIDX(1). The
%data point corresponding to this is data(:,DATAIDX(1)). The tree is split
%at a particular node according to the discriminating dimensions given in
%DISC. For example, LOSON(1) gives the index of the subtree with nodes
%whose dimension DISC(1) is less than data(DISC(1),DATAIDX(1)). That is,
%the node at the root of the subtree is DATAIDX(LOSON(1)). On the other
%hand, HISON(1) points to the node at the root of the subtree of nodes
%whose discriminating index is greater than or equal to the one at the
%split node.
%
%Modification of the CPPData member of this class or of the error checking
%code in the members can potentially lead to Matlab crashing as CPPData is
%a pointer to data in the C++ implementation.
%
%REFERENCES:
%[1] J. L. Bentley, "Multidimensional binary search trees used for
%    associative searching," Communications of the ACM, vol. 18, no. 9, pp.
%    509-517, Sep. 1975.
%[2] J. H. Friedman, J. L. Bentley, and R. A. Finkel, "An algorithm for
%    finding best matches in logarithmic expected time," ACM Transactions
%    on Mathematical Software, vol. 3, no. 3, pp. 209-226, Sep. 1977.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

properties(Access=public)
   data%A matrix of the data points.
end

properties(Access=private)
   LOSON%Lists the index of the next lower node.
   HISON%Lists the index of the next higher node.
   DATAIDX%The index of the data at a node.
   DISC%The discriminating dimension index at the node.
   subtreeSizes%Holds the size of all of the nodes in the subtree from a given node (including the given node).
   BMin%An array of bounds of the children of each node for a given level.
   BMax
   
   CPPData%Only used if an interface to a C++ implementation exists.
end

methods
    function newTree=kdTree(k,N)
    %%KDTREE Construct a new kd tree with space to hold all of the data
    %        that is to be passed.
    %
    %INPUTS: k The dimensionality of the data.
    %        N The number of data points that will be in the tree.
    %
    %OUTPUTS: newTree A new kdTree instance with the proper amount of
    %                 space.
    %
    %Once a tree has been allocated, it can be initialized using the
    %buildTreeFromBatch method. The size of the data batch given with that
    %method should match the size of the tree allocated here.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
        
        if(exist('kdTreeCPPInt','file'))
            newTree.CPPData=kdTreeCPPInt('kdTreeCPP',k,N);
        else
            newTree.LOSON=zeros(N,1);
            newTree.HISON=zeros(N,1);
            newTree.DATAIDX=zeros(N,1);
            newTree.DISC=zeros(N,1);
            newTree.subtreeSizes=zeros(N,1);
            newTree.BMin=zeros(k,N);
            newTree.BMax=zeros(k,N);
            newTree.data=[];
        end
    end
    
    function buildTreeFromBatch(theTree,dataBatch)
    %BUILDTREEFROMBATCH Build a balanced k-d tree from a batch of data.
    %                   This will destroy any previous initialization of
    %                   the tree.
    %
    %INPUTS: theTree   The implicitly passed kdTree object.
    %        dataBatch A kXN array of points that are to be stored in the
    %                  kd tree. k is the dimensionality of the points and N
    %                  is the number of points.
    %
    %The tree is generally constructed as described in [1]. with all of the
    %nodes stored in array, but with a few changes.
    %
    %Specifically, the array subtreeSizes is added, which makes the
    %rangeCount function able to execute quickly without having to visit
    %all of the nodes on a subtree. The bounds arrays, given by BMin
    %and BMax, at a particular node bound not just the values of the
    %discriminator key at the node, but rather the values of all of the
    %dimensions of the vector. Additionally, an attempt is made to build a
    %balanced kd tree by setting the split point to the median point based
    %on the discriminant dimension at each node. However, if a lot of nodes
    %have identical values in various dimensions, then the tree can be
    %unbalanced.
    %
    %REFERENCES:
    %[1] J. L. Bentley, "Multidimensional binary search trees used for 
    %    associative searching," Communications of the ACM, vol. 18, no. 9,
    %    pp. 509-517, Sep. 1975.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
        
        if(exist('kdTreeCPPInt','file'))
            N=kdTreeCPPInt('getN',theTree.CPPData);
            k=kdTreeCPPInt('getk',theTree.CPPData);
            
            if(N~=size(dataBatch,2)||k~=size(dataBatch,1))
               error('The data batch given does not match the preallocated size'); 
            end
            
            kdTreeCPPInt('buildTreeFromBatch',theTree.CPPData,dataBatch);
            
            %This makes an extra copy of the data (since it it already
            %saved in C), but it makes it easy for folks to access the
            %points after a range query or the like.
            theTree.data=dataBatch;
        else
            k=size(dataBatch,1);
            N=size(dataBatch,2);
            
            if(N~=size(dataBatch,2)||k~=size(dataBatch,1))
               error('The data batch given does not match the preallocated size'); 
            end
            
            theTree.data=dataBatch;

            %Allocate space for k-subarrays by which the data will be split.
            sortArray=dataBatch;
            idx=1:N;

            %The bounds initially encompass everything.
            theTree.BMin=inf(k,N);
            theTree.BMax=-inf(k,N);

            theTree.treeGrow(sortArray,idx,1,1);
        end
    end
    
    function retSet=rangeQuery(theTree,rectMin,rectMax)
    %%RANGEQUERY Perform one or more orthogonal range queries for the given
    %            bounds specified in rectMin and rectMax.
    %
    %INPUT: theTree The implicitly passed kdTree object.
    %       rectMin For k-dimensional data, this is a kXm vector of the
    %               minimum bounds of the m hyperrectangles in which one is
    %               to search for data points.
    %       rectMax The maximum bounds of the m hyperrectangles in which
    %               one is to search for data points.
    %
    %OUTPUTS: retSet An instance of the ClusterSet class containing all of
    %                the range query results. If no points fall into a
    %                hyperrectangle, then the corresponding clusterSize
    %                element in the ClusterSet will be zero.
    %
    %The orthogonal range query algorithm is essentially that described in
    %[1].
    %
    %EXAMPLE:
    %Here we demonstrate that the range query provides the same data as one
    %obtains via a brute-force query:
    % numPoints=300;
    % points=20*rand(3,numPoints)-10;
    % 
    % theTree=kdTree(3,numPoints);
    % theTree.buildTreeFromBatch(points);
    % 
    % rectMin=[-8;-8;-8];
    % rectMax=[-1;-1;-1];
    % 
    % %Do an orthogonal range query and extract the points.
    % retSet=theTree.rangeQuery(rectMin,rectMax);
    % boxedPoints=theTree.data(:,retSet.clusterEls);
    % 
    % %Do the range query by brute force comparison.
    % sel=(points(1,:)>rectMin(1))&(points(1,:)<rectMax(1));
    % sel=sel&(points(2,:)>rectMin(2))&(points(2,:)<rectMax(2));
    % sel=sel&(points(3,:)>rectMin(3))&(points(3,:)<rectMax(3));
    % selPoints=points(:,sel);
    % 
    % %We will now sort the points obtained by each query in order of
    % %increasing x value (we do not expect two points generated by the
    % %rand function to have the same x value; thus, this should uniquely
    % %sort the points).
    % [~,idx]=sort(boxedPoints(1,:),'ascend');
    % boxedPoints=boxedPoints(:,idx);
    % [~,idx]=sort(selPoints(1,:),'ascend');
    % selPoints=selPoints(:,idx);
    % %This should be true (and also not ahve an error due to a differing
    % %number of points used).
    % all(boxedPoints(:)==selPoints(:))
    % 
    %REFERENCES:
    %J. L. Bentley, "Multidimensional binary search trees used for 
    %associative searching," Communications of the ACM, vol. 18, no. 9, pp.
    %509-517, Sep. 1975.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.

        m1=size(rectMin,2);
        m2=size(rectMax,2);
        
        if(m1~=m2)
            error('Inconsistent numbers of rectangles given.');
        end
        
        if(exist('kdTreeCPPInt','file'))
            N=kdTreeCPPInt('getN',theTree.CPPData);
            k=kdTreeCPPInt('getk',theTree.CPPData);
            
            if(k~=size(rectMin,1)||k~=size(rectMax,1))
               error('The coordinates are not the appropriate sizes.'); 
            end

            if(N==0)
               retSet=ClusterSet([],0,0);
               return;
            end
            
            retSet=kdTreeCPPInt('rangeQuery',theTree.CPPData,rectMin,rectMax);
            %The +1 converts C indicies to Matlab indicies
            retSet.clusterEls=retSet.clusterEls+1;
        else
            %Create a new ClusterSet with m1 empty clusters.
            retSet=ClusterSet([],zeros(m1,1),zeros(m1,1));
            
            %Add the clusters
            for curRect=1:m1
                idxRange=theTree.rangeQueryRecur(1,rectMin(:,curRect),rectMax(:,curRect));
                retSet.clusterSizes(curRect)=size(idxRange,1);
                retSet.clusterEls=[retSet.clusterEls;idxRange];
                if(curRect<m1)
                    retSet.offsetArray(curRect+1)=retSet.offsetArray(curRect)+retSet.clusterSizes(curRect);
                end
            end
        end
    end
    
    function numInRange=rangeCount(theTree,rectMin,rectMax)
    %%RANGECOUNT Count the number of elements in one or more orthogonal
    %            range queries without actually returning the elements.
    %
    %INPUT: theTree  The implicitly passed kdTree object.
    %       rectMin  For k-dimensional data, this is a kXm vector of the
    %                minimum bounds of the m hyperrectangles in which one
    %                is to search for data points.
    %       rectMax The maximum bounds of the m hyperrectangles in which
    %               one is to search for data points.
    %
    %OUTPUTS: numInRange The NX1 vvector of the integer number of entries
    %                    in the kd-tree that fall into each hyperrectangle.
    %
    %This is essentially the same as the rangeQuery method, except entire
    %subtrees can be skipped, since the number of nodes in a subtree is
    %already known due to the subtreeSizes array.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        m1=size(rectMin,2);
        m2=size(rectMax,2);
        
        if(m1~=m2)
            error('Inconsistent numbers of rectangles given.');
        end
        
        if(exist('kdTreeCPPInt','file'))
            N=kdTreeCPPInt('getN',theTree.CPPData);
            k=kdTreeCPPInt('getk',theTree.CPPData);
            
            if(k~=size(rectMin,1)||k~=size(rectMax,1))
               error('The coordinates are not the appropriate sizes.'); 
            end
            
            if(N==0)
               numInRange=0;
               return;
            end
            
            numInRange=kdTreeCPPInt('rangeCount',theTree.CPPData,rectMin,rectMax);
        else
            numInRange=zeros(m1,1);
            
            for curRect=1:m1
                numInRange(curRect)=theTree.rangeCountRecur(1,rectMin(:,curRect),rectMax(:,curRect));
            end
        end
    end
    
    function [idxRange, distSquared]=findmBestNN(theTree,point,m)
    %%FINDMBESTNN Return the indices of the k-best nearest neighbors
    %             (according to the squared l2 norm) of the given point and
    %             the squared distances of the nearest neighbors to the
    %             point.
    %
    %INPUTS: theTree The implicitly passed kdTree object.
    %          point A kXn matrix of n points whose m nearest neighbors are
    %                desired.
    %              m The number of nearest neighbors to find for each
    %                point. If m > the number of elements in the k-d tree,
    %                then an error is raised.
    %
    %OUTPUTS: idxRange A kX1 vector such that theTree.data(:,idxRange(k))
    %                  is the k-best match.
    %      distSquared A kX1 vector of the squared Euclidean distance from
    %                  points found in idxRange to the given point.
    %
    %This is essentially the algorithm described in [1].
    %
    %REFERENCES:
    %[1] J. H. Friedman, J. L. Bentley, and R. A. Finkel, "An algorithm for
    %    finding best matches in logarithmic expected time," ACM
    %    Transactions on Mathematical Software, vol. 3, no. 3, pp. 209-226,
    %    Sep. 1977.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        if(nargin<3)
            m=1;
        end

        if(exist('kdTreeCPPInt','file'))
            N=kdTreeCPPInt('getN',theTree.CPPData);
            k=kdTreeCPPInt('getk',theTree.CPPData);
            kPoint=size(point,1);
            
            if(m>N)
                error('More neighbors requested than there are elements in the tree.');
            end
            
            if(k~=kPoint)
                error('The points have the wrong dimensionality.');
            end
            
            [idxRange, distSquared]=kdTreeCPPInt('findmBestNN',theTree.CPPData,point,m);
            idxRange=idxRange+1;%Convert C indicies to Matlab indicies.
        else
            N=size(theTree.data,2);
            k=size(theTree.data,1);
            numPoints=size(point,2);
            kPoint=size(point,1);
            
            if(m>N)
                error('More neighbors requested than there are elements in the tree.');
            end
            
            if(k~=kPoint)
                error('The points have the wrong dimensionality.');
            end

            idxRange=zeros(m,numPoints);
            distSquared=zeros(m,numPoints);
            for curPoint=1:numPoints
                %Allocate a priority queue with enough space for the found
                %results. The value at the start of the queue is the one
                %with the highest cost.
                mBestQueue=BinaryHeap(m);
                %Start the recursion to fill the queue with the k-best values.
                theTree.mBestRecur(1,mBestQueue,point(:,curPoint),m);

                %Extract the k-best values from the queue to return.
                mBestQueue.heapSize;

                for curFound=m:-1:1
                    topPair=mBestQueue.deleteTop();
                    distSquared(curFound,curPoint)=topPair.key;
                    idxRange(curFound,curPoint)=theTree.DATAIDX(topPair.value);
                end
                mBestQueue.delete;
            end
        end
    end
    
    function display(theTree)
    %%DISPLAY Display information about the tree, including whether the C++
    %         implementation (wrapped by a Matlab class) or the Matlab-only
    %         implementation is being used.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        if(exist('kdTreeCPPInt','file'))
            display('A kdTree as implemented via a mex wrapper around a C++ class.')
        else
            display('A kdTree as implemented as a Matlab class.')
        end
        theTree.disp();
    end
    
    function copyDataFromCPP(theTree)
    %%COPYDATAFROMCPP Copy the data that makes up the structure of the tree
    %                 out of C++ into the variables in the Matlab kdTree
    %                 class. This can be useful when debugging the C++
    %                 implementation. This function just raises an error
    %                 when executed when using the Matlab implementation.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
        
        if(exist('kdTreeCPPInt','file'))
            [theTree.LOSON,theTree.HISON,theTree.DATAIDX,theTree.DISC,theTree.subtreeSizes,theTree.BMin,theTree.BMax,theTree.data]=kdTreeCPPInt('getAllData',theTree.CPPData);
        else
            error('This kdTree class instance does not stem from a C++ class');
        end
    end
    
    function delete(theTree)
    %%DELETE The destructor method. This method is used when the kd tree is
    %        implemented as a C++ class. This method prevents a memory
    %        leak.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
        
        if(exist('kdTreeCPPInt','file'))
            kdTreeCPPInt('~kdTreeCPP',theTree.CPPData);
        end
    end
end

methods(Access=private)
    function nextFreeNode=treeGrow(theTree,sortArray,idx,level,curNode)
    %TREEGROW The recursion function for building a kd tree from a batch of
    %         data.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 
        k=size(theTree.data,1);
        NSubTree=length(idx);
        
        theTree.subtreeSizes(curNode)=NSubTree;
        %Record the discriminating index at this level.
        theTree.DISC(curNode)=level;
        
        nextFreeNode=curNode+1;
        %First, check whether this is a leaf node. If so, then it has no
        %children and the index of the discriminator is just idx.
        if(NSubTree==1)
            theTree.LOSON(curNode)=-1;
            theTree.HISON(curNode)=-1;
            theTree.DATAIDX(curNode)=idx;
            
            theTree.BMin(:,curNode)=sortArray(:,1);
            theTree.BMax(:,curNode)=sortArray(:,1);
            return;
        end
        
        %Sort the elements of sortArray and idx according to the indices of
        %the elements in the current level.
        [~,idxNew]=sort(sortArray(level,:),'ascend');
        sortArray=sortArray(:,idxNew);
        idx=idx(idxNew);
        
        %Find the first occurance of the median element.
        midIdx=findFirstMax(sortArray(level,1:floor((NSubTree+1)/2)));
        
        %The median (split) element.
        theTree.DATAIDX(curNode)=idx(midIdx);
                
        %The next level's index.
        nextLevel=mod(level,k)+1;
        
        %If duplicate values are present, there might be nothing before the
        %current node and LOSON will be empty.
        if(midIdx==1)
            theTree.LOSON(curNode)=-1;%There is no splitting.
        else
            theTree.LOSON(curNode)=nextFreeNode;
            nextNode=theTree.treeGrow(sortArray(:,(1:(midIdx-1))),idx(1:(midIdx-1)),nextLevel,nextFreeNode);
            
            %Record the minimum and maximum values.
            theTree.BMin(:,curNode)=theTree.BMin(:,nextFreeNode);
            theTree.BMax(:,curNode)=theTree.BMax(:,nextFreeNode);
            nextFreeNode=nextNode;
        end
        
        %The HISON node will never be empty, except at a leaf node, since
        %the median is always taken using the floor function.
        theTree.HISON(curNode)=nextFreeNode;
        nextNode=theTree.treeGrow(sortArray(:,(midIdx+1):end),idx((midIdx+1):end),nextLevel,nextFreeNode);
        
        %Record the minimum and maximum values due to the child node and
        %due to the contribution of the current node.
        theTree.BMin(:,curNode)=min([theTree.BMin(:,curNode),theTree.BMin(:,nextFreeNode),sortArray(:,midIdx)],[],2);
        theTree.BMax(:,curNode)=max([theTree.BMax(:,curNode),theTree.BMax(:,nextFreeNode),sortArray(:,midIdx)],[],2);
        nextFreeNode=nextNode;
    end
        
    function idxRange=rangeQueryRecur(theTree,curNode,rectMin,rectMax)
    %RANGEQUERYRECUR The recursion for performing an orthogonal range
    %                query.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        idxRange=[];
                
        %If curNode and all of its children are in the box, then return the
        %tree and all of its children.
        if(rectContainedInRect(theTree.BMin(:,curNode),theTree.BMax(:,curNode),rectMin,rectMax))
            idxRange=returnSubtreeIdx(theTree,curNode);
            return;
        end
        
        %Otherwise, return the point, if it is in the range, and all of the
        %points from subtrees that overlap.
        P=theTree.data(:,theTree.DATAIDX(curNode));
        if(inHyperrect(P,rectMin,rectMax))
            idxRange=[idxRange;theTree.DATAIDX(curNode)];
        end
        
        if(theTree.HISON(curNode)~=-1)
            childNode=theTree.HISON(curNode);
            %If points might be in the child nodes.
            if(rectsIntersect(theTree.BMin(:,childNode),theTree.BMax(:,childNode),rectMin,rectMax))
                idxRange=[idxRange;theTree.rangeQueryRecur(childNode,rectMin,rectMax)];
            end
        end
        
        if(theTree.LOSON(curNode)~=-1)
            childNode=theTree.LOSON(curNode);
            %If points might be in the child nodes.
            if(rectsIntersect(theTree.BMin(:,childNode),theTree.BMax(:,childNode),rectMin,rectMax))
                idxRange=[idxRange;theTree.rangeQueryRecur(childNode,rectMin,rectMax)];
            end
        end
    end
    
    function numInRange=rangeCountRecur(theTree,curNode,rectMin,rectMax)
    %RANGECOUNTRECUR The recursion for counting the number of results that
    %                would be obtained by performing an orthogonal range
    %                query.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
    
        numInRange=0;
                
        %If curNode and all of its children are in the box, then return the
        %number of items in the whole subtree.
        if(rectContainedInRect(theTree.BMin(:,curNode),theTree.BMax(:,curNode),rectMin,rectMax))
            numInRange=numInRange+theTree.subtreeSizes(curNode);
            return;
        end
        
        %Otherwise, return one, if the point is in the range, and all of 
        %the points from subtrees that overlap.
        P=theTree.data(:,theTree.DATAIDX(curNode));
        if(inHyperrect(P,rectMin,rectMax))
            numInRange=numInRange+1;
        end
        
        if(theTree.HISON(curNode)~=-1)
            childNode=theTree.HISON(curNode);
            %If points might be in the child nodes.
            if(rectsIntersect(theTree.BMin(:,childNode),theTree.BMax(:,childNode),rectMin,rectMax))
                numInRange=numInRange+theTree.rangeCountRecur(childNode,rectMin,rectMax);
            end
        end
        
        if(theTree.LOSON(curNode)~=-1)
            childNode=theTree.LOSON(curNode);
            %If points might be in the child nodes.
            if(rectsIntersect(theTree.BMin(:,childNode),theTree.BMax(:,childNode),rectMin,rectMax))
                numInRange=numInRange+theTree.rangeCountRecur(childNode,rectMin,rectMax);
            end
        end
    end
    
    
    function idxRange=returnSubtreeIdx(theTree,nodeIdx)
    %%RETURNSUBTREEIDX  Returns the indices of all of the data points in
    %                   the theTree.data that are in the subtree rooted at
    %                   nodeIndex.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
        idxRange=zeros(theTree.subtreeSizes(nodeIdx),1);
        
        %The first node in the ranges is going to be 
        idxRange(1)=theTree.DATAIDX(nodeIdx);
        minIdx=2;
        if(theTree.HISON(nodeIdx)~=-1)%If there are children to the right
            numSub=theTree.subtreeSizes(theTree.HISON(nodeIdx));
            idxRange(minIdx:(minIdx+numSub-1),1)=theTree.returnSubtreeIdx(theTree.HISON(nodeIdx));
            minIdx=minIdx+numSub;
        end
        
        if(theTree.LOSON(nodeIdx)~=-1)%If there are children to the left.
            numSub=theTree.subtreeSizes(theTree.LOSON(nodeIdx));
            idxRange(minIdx:(minIdx+numSub-1),1)=theTree.returnSubtreeIdx(theTree.LOSON(nodeIdx));
        end
    end
    
    function mBestRecur(theTree,curNodeIdx,mBestQueue,point,m)
    %MBESTRECUR The recursion function for finding the m-nearest neighbor
    %           points of a given point.
    %
    %December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
        
        %First, go down the path on the nearest side of the splitting
        %dimension from this point.
        splitPoint=theTree.data(:,theTree.DATAIDX(curNodeIdx));
        splitDim=theTree.DISC(curNodeIdx);
        
        if(point(splitDim)<splitPoint(splitDim))
            lIdx=theTree.LOSON(curNodeIdx);
            if(lIdx~=-1)
                theTree.mBestRecur(lIdx,mBestQueue,point,m);
            end
            
            farNode=theTree.HISON(curNodeIdx);
        else
            hIdx=theTree.HISON(curNodeIdx);
            if(hIdx~=-1)
                theTree.mBestRecur(hIdx,mBestQueue,point,m);
            end
            
            farNode=theTree.LOSON(curNodeIdx);
        end
        
        %Next, visit this node.
        cost=dist(point,splitPoint);
        
        if(mBestQueue.heapSize==0)
            %Since no nodes have been found to this point, the best node is
            %the current one found thus far. This path will only be visited
            %at a leaf of the tree, so we can just return.
            mBestQueue.insert(cost,curNodeIdx);
            return
        else
            %The point is only added to the queue if it is lower than the
            %maximum cost point found thus far or if there are fewer than k
            %points already in the queue.
            keyValPair=mBestQueue.getTop();
            maxDist=keyValPair.key;
            if(maxDist>cost||mBestQueue.heapSize<m)
                if(mBestQueue.heapSize==m)
                    mBestQueue.deleteTop;
                end
                mBestQueue.insert(cost,curNodeIdx);
                maxDist=cost;
            end
            
            %Now, see if it is necessary to visit the other branch of the
            %tree. That is only the case if the bounding box intersects
            %with a ball centered at the point to find whose squared radius
            %is equal to maxDist or if there are fewer than k things in the
            %queue.
            if(farNode~=-1)
                if(mBestQueue.heapSize<m||boundsIntersectBall(point,maxDist,theTree.BMin(:,farNode),theTree.BMax(:,farNode)))
                    theTree.mBestRecur(farNode,mBestQueue,point,m);
                end
            end
        end       
    end
    
end
end

function val=dist(a,b)
    diff=a-b;
    val=dot(diff,diff);
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
