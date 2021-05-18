classdef metricTree < handle
%%METRICTREE A metric tree class for performing searches in a radius
%            around a given point. If a C++ class interface for the metric
%            tree has been compiled, then its functions will be called in
%            place of the Matlab routines, since the C++ implementation
%            can be significantly faster.
%
%The metric tree data structure is implemented as described in [1] and [2].
%
%This implementation builds the tree from a batch of data all at once. As
%long as there are not numerous points equidistance from a given point, the
%tree will be balanced. The Euclidean distance is used as the metric,
%though the code could be changed to replace it with almost any distance 
%metric that satisfies the triangle inequality.
%
%Note that unlike certain metric tree implementations, not all of the nodes
%are held in the leaves.
%
%Modification of the CPPData member of this class or of the error checking
%code in the members can potentially lead to Matlab crashing as CPPData is
%a pointer to data in the C++ implementation.
%
%REFERENCES:
%[1] J. K. Uhlmann. (1991, Nov.) Implementing metric trees to satisfy
%    general proximity/similarity queries. [Online]. Available:
%    http://people.cs.missouri.edu/ uhlmannj/ImplementGH.pdf
%[2] J. K. Uhlmann, "Satisfying general proximity/similarity queries with
%    metric trees," Information Processing Letters, vol. 40, no. 4, pp.
%    175-179, 25 Nov. 1991.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    properties
        DATAIDX%The index of the data at a node.
        innerChild
        outerChild
        innerRadii
        outerRadii
        data%A matrix of the data points.
        
        CPPData%Only used if an interface to a C++ implementation exists.
    end
   
    methods
        function newTree=metricTree(k,N)
        %%METRICTREE Construct a new metric tree with space to hold a given
        %            number of nodes.
        %
        %INPUTS: k The dimensionality of the nodes that will be placed in
        %          the tree.
        %        N The number of nodes that the tree will hold.
        %
        %OUTPUTS: newTree A new metricTree instance with the proper amount
        %                 of space.
        %
        %Once a tree has been allocated, it can be initialized using the
        %buildTreeFromBatch method. The size of the data batch given with
        %that method should match the size of the tree allocated here.
            if(exist('metricTreeCPPInt','file'))
                newTree.CPPData=metricTreeCPPInt('metricTreeCPP',k,N);
            else
                newTree.DATAIDX=zeros(N,1);
                newTree.innerChild=zeros(N,1);
                newTree.outerChild=zeros(N,1);
                newTree.innerRadii=zeros(N,1);
                newTree.outerRadii=zeros(N,1);
                newTree.data=[];
            end
        end
        
        function buildTreeFromBatch(theTree,dataBatch)
        %BUILDTREEFROMBATCH Build a metric tree from a batch of data.
        %
        %INPUTS: theTree The implicitly passed metricTree object.
        %      dataBatch A kXN array of points that are to be stored in
        %                the metric tree. k is the dimensionality of the
        %                points and N is the number of points.
        %
        %The tree is constructed roughly as described in [1].
        %
        %REFERENCES:
        %[1] J. K. Uhlmann. (1991, Nov.) Implementing metric trees to
        %    satisfy general proximity/similarity queries. [Online].
        %    Available: http://people.cs.missouri.edu/ uhlmannj/ImplementGH.pdf
           
            N=size(dataBatch,2);
            if(exist('metricTreeCPPInt','file'))
                N=metricTreeCPPInt('getN',theTree.CPPData);
                k=metricTreeCPPInt('getk',theTree.CPPData);
            
                if(N~=size(dataBatch,2)||k~=size(dataBatch,1))
                   error('The data batch given does not match the preallocated size'); 
                end
                
                metricTreeCPPInt('buildTreeFromBatch',theTree.CPPData,dataBatch);
            else
                theTree.data=dataBatch;

                %Preallocate an adjacency matrix that will be used in the
                %recursion.
                adjMat=zeros(N,1);
                idx=(1:N)';
                theTree.treeGrow(adjMat,dataBatch,idx,1);
            end
        end
       
        function [retSet,distSet]=searchRadius(theTree,point,radius)
        %%SEARCHRADIUS Return all of the points in the metric tree that are
        %              located within a given search radius of within a
        %              given radius about a given point.
        %
        %INPUTS: theTree The implicitly passed metricTree object.
        %          point A kXm matrix of k-dimensional points about which
        %                range searches will be performed.
        %         radius An mX1 vector of the search radii about each
        %                point.
        %
        %OUTPUTS: retSet An instance of the ClusterSet class containing all
        %                of the items within the given radius for each
        %                point. If no points in the tree fall into a given
        %                region, then the corresponding clusterSize element
        %                in the ClusterSet will be zero.
        %        distSet An instance of the clusterSet class containing all
        %                of the distances associated with the found items
        %                in retSet

            m1=size(point,2);
            m2=size(radius,1);

            if(m1~=m2)
                error('Inconsistent numbers of points and radii given.');
            end
            
            if(exist('metricTreeCPPInt','file'))
                [retSet,distSet]=metricTreeCPPInt('searchRadius',theTree.CPPData,point,radius);
                %Convert indices to Matlab indices
                retSet.clusterEls=retSet.clusterEls+1;
            else
                %Create the new ClusterSet instances with m1 empty clusters.
                retSet=ClusterSet([],zeros(m1,1),zeros(m1,1));
                distSet=ClusterSet([],zeros(m1,1),zeros(m1,1));

                %Add the clusters
                for curPoint=1:m1
                    [idxRange,distList]=theTree.searchRadRecur(point(:,curPoint),radius(curPoint),1);

                    retSet.clusterSizes(curPoint)=size(idxRange,1);
                    retSet.clusterEls=[retSet.clusterEls;idxRange];
                    distSet.clusterSizes(curPoint)=size(idxRange,1);
                    distSet.clusterEls=[distSet.clusterEls;distList];
                    if(curPoint<m1)
                        retSet.offsetArray(curPoint+1)=retSet.offsetArray(curPoint)+retSet.clusterSizes(curPoint);
                        distSet.offsetArray(curPoint+1)=distSet.offsetArray(curPoint)+distSet.clusterSizes(curPoint);
                    end
                end
            end
        end
        
        function display(theTree)
        %%DISPLAY Display information about the tree, including whether the
        %         C++ implementation (wrapped by a Matlab class) or the
        %         Matlab-only implementation is being used.
    
            if(exist('metricTreeCPPInt','file'))
                display('A metric tree as implemented via a mex wrapper around a C++ class.')
            else
                display('A metric tree as implemented as a Matlab class.')
            end
            theTree.disp();
        end
        
        function copyDataFromCPP(theTree)
    %%COPYDATAFROMCPP Copy the data that makes up the structure of the tree
    %                 out of C++ into the variables in the Matlab kdTree
    %                 class. This can be useful when debugging the C++
    %                 implementation. This function just raises an error
    %                 when executed when using the Matlab implementation.
        
            if(exist('metricTreeCPPInt','file'))
                [theTree.DATAIDX,theTree.innerChild,theTree.outerChild,theTree.innerRadii,theTree.outerRadii,theTree.data]=metricTreeCPPInt('getAllData',theTree.CPPData);
            else
                error('This metric tree class instance does not stem from a C++ class');
            end
        end
    
        function delete(theTree)
        %%DELETE The destructor method. This method is used when the metric
        %        tree is implemented as a C++ class. This method prevents a
        %        memory leak.

            if(exist('metricTreeCPPInt','file'))
                metricTreeCPPInt('~metricTreeCPP',theTree.CPPData);
            end
        end
    end
    
    methods(Access=private)
        function [idxRange,distList]=searchRadRecur(theTree,point,r2,curNode)
        %%SEARCHRADRECUR A recursion function for performing a search of a
        %                particular radius about a point.
        
            distCur=dist(point,theTree.data(:,theTree.DATAIDX(curNode)));
            if(distCur<=r2)
                idxRange=theTree.DATAIDX(curNode);
                distList=distCur;
            else
                idxRange=[];
                distList=[];
            end
            
            if(distCur+r2>=theTree.outerRadii(curNode)&&theTree.outerChild(curNode)~=-1)
                [idxRest,distSqRest]=theTree.searchRadRecur(point,r2,theTree.outerChild(curNode));
                idxRange=[idxRange;idxRest];
                distList=[distList;distSqRest];
            end
            
            if(distCur-r2<=theTree.innerRadii(curNode)&&theTree.innerChild(curNode)~=-1)
                [idxRest,distSqRest]=theTree.searchRadRecur(point,r2,theTree.innerChild(curNode));
                idxRange=[idxRange;idxRest];
                distList=[distList;distSqRest];
            end
        end
        
        function nextFreeNode=treeGrow(theTree,adjMat,dataBatch,idx,curNode)
        %%TREEGROW A recursion function for creating a new metric tree.
            
            NSubTree=size(dataBatch,2);
            
            %The first index in idx is added to the tree.
            theTree.DATAIDX(curNode)=idx(1);
            
            %If this is a leaf node.
            if(NSubTree==1)
                theTree.innerChild(curNode)=-1;
                theTree.outerChild(curNode)=-1;
                theTree.innerRadii(curNode)=-Inf;
                theTree.outerRadii(curNode)=Inf;
                nextFreeNode=curNode+1;
                return
            end

            %Fill the NCurX1 subvector of adjMat with all pairwise
            %distances from the current point to the other points in idx.
            NSubTree=NSubTree-1;
            for curPoint=1:NSubTree
                adjMat(curPoint)=dist(dataBatch(:,1),dataBatch(:,curPoint+1));
            end
            
            %Next, sort the points
            [adjMat(1:NSubTree), sortIdx]=sort(adjMat(1:NSubTree));
            
            idx=idx(2:end);
            idx=idx(sortIdx);
            
            dataBatch=dataBatch(:,2:end);
            dataBatch=dataBatch(:,sortIdx);
            
            %Find the first occurance of the median distance.
            midIdx=findFirstMax(adjMat(1:ceil((NSubTree+1)/2)));
            
            %The inner and outer radii of how things are split in terms of
            %distances from the given point.
            theTree.outerRadii(curNode)=adjMat(midIdx);
            if(midIdx>1)
                theTree.innerRadii(curNode)=adjMat(midIdx-1);
            else
                theTree.innerRadii(curNode)=-Inf;
            end
            
            %Continue the recursion.
            theTree.outerChild(curNode)=curNode+1;
            nextFreeNode=theTree.treeGrow(adjMat,dataBatch(:,midIdx:end),idx(midIdx:end),curNode+1);
            
            %If a full partitioning of the nodes took place
            if(midIdx>1)
                theTree.innerChild(curNode)=nextFreeNode;
                nextFreeNode=theTree.treeGrow(adjMat,dataBatch(:,1:(midIdx-1)),idx(1:(midIdx-1)),nextFreeNode);
            else
                theTree.innerChild(curNode)=-1;
            end
        end
    end
end

function val= dist(a,b)
%%DIST Evaluate the Euclidean distance between two vectors.

    %The square root must be used to obey the triangle inequality.
    diff=a-b;
    val=sqrt(diff'*diff);
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
