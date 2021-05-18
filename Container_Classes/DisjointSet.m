classdef DisjointSet < handle
%%DISJOINTSET A disjoint set class that can be used for partitioning a
%             general set. It can be used to cluster targets given
%             measurements. If one also wishes to keep track of which
%             measurements are associated with which target clusters, then
%             the subclass DisjointSetM should be used.
%
%Disjoint sets are described in Chapter 8 of [1], in Chapter 21 of [2], and
%Chapter 12.4 of [3].
%
%The main point of this class is to bring together cluster data so that a
%ClusterSet class can be created. A ClusterSet class allows the clusters to
%be indexed by a cluster number and the elements of the clusters to be
%subsequently indexed with another number. The algorithms here use path
%compression in the find algorithm and a union by rank algorithm to process
%things quickly.
%
%REFERENCES:
%[1] M. A. Weiss, Data Structures and Algorithm Analysis in C++, 2nd ed.
%    Reading, MA: Addison-Wesley, 1999.
%[2] T. H. Cormen, C. E. Leiserson, R. L. Rivest, and C. Stein,
%    Introduction to Algorithms, 2nd ed. Cambridge, MA: The MIT Press,
%    2001.
%[3] B. R. Preiss, Data Structures and Algorithms with Object-Oriented
%    Design Patterns in C++. New York, NY: John Wiley & Sons, Inc., 1999.
%
%November 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        numberOfSets
        setArray
        numInTree
    end
    methods
        function newSet=DisjointSet(setSize)
            %The constructor method
            newSet.setArray=-1*ones(setSize,1);
            newSet.numInTree=ones(setSize,1);
            newSet.numberOfSets=setSize;
        end
        
        function unionFromList(DSObj,uList)
        %%UNIONFROMLIST Given a list of targets in the same cluster, such
        %               as a list of targets that gate with a common
        %               measurement, cluster them together.
        %
        %INPUTS: DSObj The implicitly passed calling object.
        %        uList A one-dimensional list of indices of targets to
        %              merge into one cluster
        
            numInList=length(uList);
            if(numInList==0)
                return;
            end
            
            rootIdx1=DSObj.find(uList(1));
            for curTar=2:numInList
                rootIdx2=DSObj.find(uList(curTar));
                rootIdx1=DSObj.unionRoots(rootIdx1,rootIdx2);
            end
        end
        
        function unionFromBinMat(DSObj,binMat)
        %%UNIONFROMBINMAT Given a binary matrix specifying which 
        %                 measurements gate with which targets, cluster
        %                 together all targets that gate with each other.
        %
        %INPUTS: DSObj The implicitly passed calling object.
        %        binMat A numTarXnumMeas binary matrix where numTar is
        %               the same as the length of setArray in DSObj and
        %               numMeas is an arbitrary number. If an empty matrix
        %               is passed, nothing is done.
        %
        %This algorithm takes the union of all elements of setArray in each
        %column corresponding to ones in the column. Thus, if this is a
        %gating matrix for targets and measurements, this will create all
        %of the target clusters.
        %
            if(isempty(binMat))
                return;
            end

            numTar=size(binMat,1);
            numMeas=size(binMat,2);

            %Check for input validity
            assert(numTar==length(DSObj.setArray))
            
            for curMeas=1:numMeas
                rootIdx1=0;
                for curTar=1:numTar
                    if(binMat(curTar,curMeas)~=0)
                        if(rootIdx1==0)
                            rootIdx1=DSObj.find(curTar);
                        else
                            rootIdx2=DSObj.find(curTar);
                            rootIdx1=DSObj.unionRoots(rootIdx1,rootIdx2);
                        end
                    end
                end
            end
        end

        function newCSet=createClusterSet(DSObj)
        %%CREATECLUSTERSET Turn the DisjointSet object into a ClusterSet
        %                  so that targets in each cluster can be easily
        %                  addressed.
        %
        %OUTPUTS: newCSet A ClusterSet object representing the clustering
        %                 given by this DisjointSet object.
        
            numTar=length(DSObj.setArray);
            
            rootIdxList=zeros(DSObj.numberOfSets,1);
            curNumAdded=zeros(DSObj.numberOfSets,1);
            
            rootRankList=zeros(DSObj.numberOfSets,1);
            clusterEls=zeros(numTar,1);
            
            %First, record the indices and ranks of all of the root nodes.
            rootNodesFound=0;
            for curItem=1:length(DSObj.setArray)
                if(DSObj.setArray(curItem)<0)
                    rootNodesFound=rootNodesFound+1;
                    rootIdxList(rootNodesFound)=curItem;
                    rootRankList(rootNodesFound)=DSObj.numInTree(curItem);
                    if(rootNodesFound==DSObj.numberOfSets)
                        break;
                    end
                end
            end
            
            %Next, sort the root nodes by their index, so that children can
            %be efficiently assigned to them.
            [rootIdxList,idx] = sort(rootIdxList,1,'ascend');
            rootRankList=rootRankList(idx);
            
            %rootRankList currently holds all of the ranks of the root
            %indices in rootIdxList for the clusters. The elements of the
            %clusters are going to be put in one big array, clusterEls.
            %We are going to  create another array, offsetArray from
            %rootRankList such that each element is the offset needed to
            %access the first element of a given cluster in clusterEls. The
            %first element of offsetArray consequently is zero. The others
            %are a cumulative sum of the previous elements in rootRankList.
            offsetArray=zeros(DSObj.numberOfSets,1);
            offsetArray(2:end)=cumsum(rootRankList(1:(end-1)));
            
            %Now, we shall scan though all of the elements in
            %DSObj.setArray and add then to the appropriate cluster in
            %clusterEls depending on their rank.
            for curItem=1:numTar
                %Find the root node for the given node.
                rootNode=DSObj.find(curItem);
                %Find the index of the cluster corresponding to that root
                %node.
                [~, clusterIdx]=binSearch(rootIdxList,rootNode);
                %Add the node to the cluster.
                clusterEls(1+offsetArray(clusterIdx)+curNumAdded(clusterIdx))=curItem;
                curNumAdded(clusterIdx)=curNumAdded(clusterIdx)+1;
            end
            
            newCSet=ClusterSet(clusterEls,rootRankList,offsetArray);
        end

        function mergedRootIdx=unionRoots(DSObj,rootIdx1,rootIdx2)
        %%UNIONROOTS Merge two sets given their root nodes. If one wishes
        %            to merge sets given arbitrary (not necessarily root)
        %            nodes in the sets, then the function union should be
        %            used.
        %
        %INPUTS: DSObj      The implicitly passed calling object.
        %        rootIdx1   The root index of the first set to merge.
        %        rootIdx2   The root index of the second set to merge.
        %
        %OUTPUTS: mergedRootIdx The index of the root of the merged
        %                       cluster. It will either be rootIdx1 or
        %                       rootIdx2.
        %
        %This is a union-by-rank algorithm.
        %
            %If the roots are the same, then there is nothing to merge.
            if(rootIdx1==rootIdx2)
                %Do measurement bookeeping, if desired.
                mergedRootIdx=rootIdx1;
                return;
            end
            
            %If root 2 is deeper, than make root2 the new root.
            if(DSObj.setArray(rootIdx2)< DSObj.setArray(rootIdx1))
                DSObj.setArray(rootIdx1)=rootIdx2;
                DSObj.numInTree(rootIdx2)=DSObj.numInTree(rootIdx1)+DSObj.numInTree(rootIdx2);
                mergedRootIdx=rootIdx2;
            else
                %If they are the same height, then update the height of
                %rootIdx1.
                if(DSObj.setArray(rootIdx1)== DSObj.setArray(rootIdx2))
                    DSObj.setArray(rootIdx1)=DSObj.setArray(rootIdx1)-1;
                end
                DSObj.setArray(rootIdx2)=rootIdx1;
                DSObj.numInTree(rootIdx1)=DSObj.numInTree(rootIdx1)+DSObj.numInTree(rootIdx2);
                mergedRootIdx=rootIdx1;
            end

            %Merging reduces the number of sets by one.
            DSObj.numberOfSets=DSObj.numberOfSets-1;
        end

        function rootIdx=find(DSObj,idx)
        %FIND Perform a search with path compression. That is, all nodes
        %     traversed are set to point directly to the root of the tree.
        %     This means that future searches will be faster.
        %
        %INPUTS:  DSObj     The implicitly passed calling object.
        %         idx       The index of a node.
        %
        %OUTPUTS: rootIdx  The index of the root node of the cluster to
        %                  which the node with index idx belongs.
        %
        %The use of path compression would imply that the ranks associated
        %with the root nodes should be recalculated. However, it is noted
        %in Chapters 8.5 and 8.6 of 
        %M. A. Weiss, Data Structures and Algorithms for C++, 2nd ed.
        %Reading, MA: Addison-Wesley, 1999.
        %that failing to update the rank does not worsen the theoretical
        %worst-case run time.
        %
            if(DSObj.setArray(idx)<0)
                rootIdx=idx;
            else
                rootIdx=DSObj.find(DSObj.setArray(idx));
                DSObj.setArray(idx)=rootIdx;
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
