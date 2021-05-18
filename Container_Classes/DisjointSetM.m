classdef DisjointSetM < DisjointSet
%%DISJOINTSETM A disjoint set class that can be used to partition a set of
%              targets into clusters while keeping track of which
%              measurements are associated with which cluster.
%
%This subclass generalizes the DisjointSet class to keep track of which
%measurements caused the partitioning of the targets. The set size is the
%number of ordered targets and a target is identifiable by its index number
%in the set. The class is created for a fixed number of targets and
%measurements. The main point of this class is to facilitate the creation
%of ClusterSet classes for the targets that cluster together and for the
%measurements that are associated with each cluster.
%
%November 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

properties
    numMeasInTree
    measToNode
end
methods
    function newSet=DisjointSetM(numTar,numMeas)
        %The constructor method.
        newSet=newSet@DisjointSet(numTar);
        newSet.numMeasInTree=zeros(numTar,1);
        newSet.measToNode=zeros(numMeas,1);
    end

    function unionFromList(DSObj,uList,measNum)
    %%UNIONFROMLIST Merge clusters of targets from a given common
    %               measurement. NOTE: Calling unionFromList twice with
    %               the same measurement number will produce incorrect
    %               results.           
    %
    %INPUTS: DSObj  The implicitly passed calling object.
    %        uList  A one-dimensional list of indices of targets to
    %               merge into one cluster
    %       measNum The index of the measurement.
    %
    %The uList passed for a particular measNum should a list of the
    %indices of ALL of the targets that gate with the measurement.
    %
    %November 2013 David F. Crouse, Naval Research Laboratory, Washington
    %D.C.
    
        numInList=length(uList);
        if(numInList==0)
            return;
        end

        rootIdx1=DSObj.find(uList(1));
        DSObj.numMeasInTree(rootIdx1)=DSObj.numMeasInTree(rootIdx1)+1;
        DSObj.measToNode(measNum)=rootIdx1;
        for curTar=2:numInList
            rootIdx2=DSObj.find(uList(curTar));
            rootIdx1=DSObj.unionRoots(rootIdx1,rootIdx2);
        end
    end

    function unionFromBoolList(DSObj,boolList,measNum)
    %%UNIONFROMBOOLLIST Merge targets specified by a boolean array for
    %               a given measurement index. NOTE: Calling
    %               unionFromBoolList twice with the same measurement
    %               number will produce incorrect results.    
    %
    %INPUTS: DSObj  The implicitly passed calling object.
    %     boolList  A length numTarX1 or 1XnumTar boolean array
    %               indicating whether which targets are in the same
    %               cluster due to the specified measurement.
    %       measNum The index of the measurement.
    %
    %November 2013 David F. Crouse, Naval Research Laboratory, Washington
    %D.C.
    
        numTar=length(DSObj.setArray);
        %Check for input validity.
        assert(numTar==length(DSObj.setArray))

        tarIdxCur=find(boolList,1);

        if(isempty(tarIdxCur))
            %If no target gates with the specified measurement.
            return;
        end

        rootIdx1=DSObj.find(tarIdxCur);
        DSObj.numMeasInTree(rootIdx1)=DSObj.numMeasInTree(rootIdx1)+1;
        DSObj.measToNode(measNum)=rootIdx1;
        for curTarIdx=(tarIdxCur+1):numTar
            if(boolList(curTarIdx))
                rootIdx2=DSObj.find(curTarIdx);
                rootIdx1=DSObj.unionRoots(rootIdx1,rootIdx2);
            end
        end
    end

    function unionFromBinMat(DSObj,binMat)
    %UNIONFROMBINMAT Given a binary matrix specifying which 
    %                measurements gate with which targets, cluster
    %                together all targets that gate with each other.
    %
    %INPUTS: DSObj The implicitly passed calling object.
    %       binMat A numTarXnumMeas binary matrix where numTar is the
    %              same as the length of setArray in DSObj and numMeas
    %              is the number of measurements with which the
    %              DisjointSet object was created.
    %
    %The measurement number is assumed from the column index when
    %clustering. This function should not be called multiple times with
    %different binary matrices as it will produce incorrect results.
    %Rather, multiple binary matrices can be OR-ed together prior to
    %calling this method.
    %
    %November 2013 David F. Crouse, Naval Research Laboratory, Washington
    %D.C.

        numTar=size(binMat,1);
        numMeas=size(binMat,2);

        %Check for input validity
        assert(numTar==length(DSObj.setArray))
        if(~isempty(DSObj.measToNode))
            assert(numMeas==length(DSObj.measToNode))
        end

        for curMeas=1:numMeas
            rootIdx1=0;
            for curTar=1:numTar
                if(binMat(curTar,curMeas)~=0)
                    if(rootIdx1==0)
                        rootIdx1=DSObj.find(curTar);

                        %Do measurement bookeeping.
                        DSObj.measToNode(curMeas)=curTar;
                        DSObj.numMeasInTree(rootIdx1)=DSObj.numMeasInTree(rootIdx1)+1;
                    else
                        rootIdx2=DSObj.find(curTar);
                        rootIdx1=DSObj.unionRoots(rootIdx1,rootIdx2);
                    end

                    %Do measurement bookeeping.
                    DSObj.measToNode(curMeas)=rootIdx1;
                end
            end
        end
    end

    function [newCSet,newCSMeasSet,ungatedMeas]=createClusterSet(DSObj)
    %%CREATECLUSTERSET Turn the DisjointSet object into a ClusterSet so
    %                  that targets in each cluster can be easily
    %                  addressed. Also produce a corresponding
    %                  ClusterSet for the clusters of measurements that
    %                  are associated with each cluster of targets and
    %                  produce a list of which measurements do not gate
    %                  with anything.
    %
    %OUTPUTS: newCSet A ClusterSet object representing the clustering
    %                 given by this DisjointSetM object.
    %    newCSMeasSet A ClusterSet object of the measurements that go
    %                 with each cluster of targets.
    %     ungatedMeas An array of indices of measurements that do not
    %                 gate with anything.
    %
    %newCSet(c,:) is the set of targets in cluster c and
    %newCSMeasSet(c,:) is the corresponding set of measurements.
    %
    %November 2013 David F. Crouse, Naval Research Laboratory, Washington
    %D.C.

        numTar=length(DSObj.setArray);
        numMeas=length(DSObj.measToNode);

        rootIdxList=zeros(DSObj.numberOfSets,1);
        curNumAdded=zeros(DSObj.numberOfSets,1);

        %This holds the numbers of targets per cluster.
        rootRankList=zeros(DSObj.numberOfSets,1);
        %This holds the numbers of measurements per cluster.
        rootRankListMeas=zeros(DSObj.numberOfSets,1);
        clusterEls=zeros(numTar,1);
        clusterElsMeas=zeros(numMeas,1);

        %First, record the indices and ranks of all of the root nodes.
        rootNodesFound=0;
        for curItem=1:numTar
            if(DSObj.setArray(curItem)<0)
                rootNodesFound=rootNodesFound+1;
                rootIdxList(rootNodesFound)=curItem;
                rootRankList(rootNodesFound)=DSObj.numInTree(curItem);
                rootRankListMeas(rootNodesFound)=DSObj.numMeasInTree(curItem);

                if(rootNodesFound==DSObj.numberOfSets)
                    break;
                end
            end
        end

        %Next, sort the root nodes by their index, so that children can
        %be efficiently assigned to them.
        [rootIdxList,idx] = sort(rootIdxList,1,'ascend');
        rootRankList=rootRankList(idx);
        %The measurement clusters must correspond to the target
        %clusters.
        rootRankListMeas=rootRankListMeas(idx);

        offsetArray=zeros(DSObj.numberOfSets,1);
        offsetArray(2:end)=cumsum(rootRankList(1:(end-1)));
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
        %Create the cluster set for the measurements
        newCSet=ClusterSet(clusterEls,rootRankList,offsetArray);
        %Create an array of the maximum possible size to hold the
        %ungated measurements.
        ungatedMeas=zeros(numMeas,1);
        numUngatedMeas=0;

        %Fill the measurement clusters in the same manner as done
        %for the regular clusters.
        curNumAdded=zeros(DSObj.numberOfSets,1);
        offsetArrayMeas=zeros(DSObj.numberOfSets,1);
        offsetArrayMeas(2:end)=cumsum(rootRankListMeas(1:(end-1)));
        for curItem=1:numMeas
            if(DSObj.measToNode(curItem)>0)
                rootNode=DSObj.find(DSObj.measToNode(curItem));
                %Find the index of the cluster corresponding to that root
                %node.
                [~, clusterIdx]=binSearch(rootIdxList,rootNode);
                %Add the node to the cluster.
                clusterElsMeas(1+offsetArrayMeas(clusterIdx)+curNumAdded(clusterIdx))=curItem;
                curNumAdded(clusterIdx)=curNumAdded(clusterIdx)+1;
            else
                numUngatedMeas=numUngatedMeas+1;
                ungatedMeas(numUngatedMeas)=curItem;
            end
        end

        %Shorten the array to the actual number of ungated
        %measurements.
        ungatedMeas=ungatedMeas(1:numUngatedMeas);

        newCSMeasSet=ClusterSet(clusterElsMeas,rootRankListMeas,offsetArrayMeas);
    end

    function mergedRootIdx=unionRoots(DSObj,rootIdx1,rootIdx2)
    %%unionRoots Merge two sets given their root nodes, keeping
    %            track of the number of measurements associated
    %            with both sets.
    %
    %INPUTS: DSObj      The implicitly passed calling object.
    %        rootIdx1   The root index of the first set to merge.
    %        rootIdx2   The root index of the second set to merge.
    %
    %OUTPUTS: mergedRootIdx The index of the root of the merged
    %                       cluster. It will either be rootIdx1 or
    %                       rootIdx2.
    %
    %This is a union-by-rank algorithm. The merged cluster has a number
    %of associated measurements equal to the sum of the number of
    %measurements in both clusters. This method should not be called
    %with indices that do not correspond to root nodes.
    %
    %November 2013 David F. Crouse, Naval Research Laboratory, Washington
    %D.C.

        %If the roots are the same, then there is nothing to merge.
        if(rootIdx1==rootIdx2)
            %Measurement bookeeping
            DSObj.numMeasInTree(rootIdx1)=DSObj.numMeasInTree(rootIdx1);
            mergedRootIdx=rootIdx1;
            return;
        end

        %If root 2 is deeper, than make root2 the new root.
        if(DSObj.setArray(rootIdx2)< DSObj.setArray(rootIdx1))
            DSObj.setArray(rootIdx1)=rootIdx2;
            DSObj.numInTree(rootIdx2)=DSObj.numInTree(rootIdx1)+DSObj.numInTree(rootIdx2);
            mergedRootIdx=rootIdx2;

            DSObj.numMeasInTree(rootIdx2)=DSObj.numMeasInTree(rootIdx1)+DSObj.numMeasInTree(rootIdx2);
        else
            %If they are the same height, then update the height of
            %rootIdx1.
            if(DSObj.setArray(rootIdx1)== DSObj.setArray(rootIdx2))
                DSObj.setArray(rootIdx1)=DSObj.setArray(rootIdx1)-1;
            end
            DSObj.setArray(rootIdx2)=rootIdx1;
            DSObj.numInTree(rootIdx1)=DSObj.numInTree(rootIdx1)+DSObj.numInTree(rootIdx2);
            mergedRootIdx=rootIdx1;

            DSObj.numMeasInTree(rootIdx1)=DSObj.numMeasInTree(rootIdx2)+DSObj.numMeasInTree(rootIdx1);
        end

        %Merging reduces the number of sets by one.
        DSObj.numberOfSets=DSObj.numberOfSets-1;
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
