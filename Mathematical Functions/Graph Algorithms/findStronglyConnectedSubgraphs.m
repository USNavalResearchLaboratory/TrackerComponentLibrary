function [setStartList,setLengths,setIdxList,subgraphList]=findStronglyConnectedSubgraphs(adjMat)
%%FINDSTRONGLYCONNECTEDSUBGRAPHS Given an adjacency matrix of a directed
%               graph, identify all of the strongly connected subgraphs.
%               Two vertices v and w are strongly connected if one can
%               follow the graph to go from v->w as well as from w->v. For
%               example, all cycles are strongly connected. A strongly
%               connected subgraph is one where all of the vertices in the
%               subgraph are strongly connected.
%
%INPUTS: adjMat An adjacency matrix of a graph. Vertex i is connected to
%               vertex j if adjMat(i,j)>=-Inf. adjMat(i,j)=Inf means that
%               the vertices are not connected.
%
%OUTPUTS: setStartList A numSubgraphsX1 list of the starting index of the 
%               vertices for each of the stringly connected subgraphs in 
%               setIdxList. The length of setStartList equals the number of
%               strongly connected subgraphs.
%    setLengths A numSubgraphsX1 list of the number of vertices in each
%               subgraph.
%    setIdxList A list of indices of the vertices in a subgraph. The
%               vertices in the ith subgraph are given by
%               setIdxList(setStartList(i):(setStartList(i)+setLengths(i)-1))
%  subgraphList This is a length-numSubgraphsX1 cell array where
%               subgraphList{i} contains the subset of the adjacency matrix
%               containing only the vertices of the ith strongly connected
%               subgraph. The ordering of the rows and columns of the
%               subgraph corresponds to the ordering of the indices in
%               setIdxList that list the nodes in the subgraph.
%
%The algorithm of Cheriyan and Mehlhorn [1] is implemented here. The
%complexity would be O(numVertices+numEdges) if the graph were sparse.
%however, since adjacency matrix is scanned for edges, rather
%than edges being given in a separate list, the complexity is higher. 
%
%Many authors refer to the algorithm of [2], which is essentially the same
%as that of [1], except the pseudocode in the paper contains some errors
%(e.g. popping an empty stack) and will not work as written.
%
%Consider the following example:
% adjMat=[Inf, 1,     1,      Inf,    Inf,    Inf;
%         Inf, Inf,   1,      1,      Inf,    Inf;
%         Inf, Inf,   Inf,    Inf,    Inf,    Inf;
%         Inf, Inf,   1,      Inf,    1,      Inf;
%         Inf, 1,     Inf,    Inf,    Inf,    1;
%         Inf, Inf,   1,      1,      Inf,    Inf];
% [setStartList,setLengths,setIdxList]=findStronglyConnectedSubgraphs1(adjMat)
%The results are setStartList=[1;2;6], setLengths=[1;4;1], and
%setIdxList =[3;2;6;5;4;1]; This corresponds to three stringly connected
%subgraphs. One is just vertex 3, another is just vertex 1, and the last
%one is vertices 2,6,5, and 4.
%
%REFERENCES:
%[1] J. Cheriyan and K. Mehlhorn, "Algorithms for dense graphs and
%    networks on the random access computer," Algorithmica, vol. 15, no. 6,
%    pp. 521-549, Jun. 1996.
%[2] H. N. Gabow, "Path-based depth-first search for strong and biconnected
%    components," Information Processing Letters, vol. 74, no. 3-4,
%    pp.107-114, May 2000.
%
%July 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numVertices=size(adjMat,1);

numInUnfinished=0;
unfinished=zeros(numVertices,1);%Array used as a stack.

numInRoots=0;
roots=zeros(numVertices,1);%Array used as a stack

inUnfinished=false(numVertices,1);
dfsNum=zeros(numVertices,1);
reached=false(numVertices,1);
count1=0;

%These extra arrays are so that we can keep track of the sets that were
%found; they are not explicitly described in the paper.
setStartList=zeros(numVertices,1);
setIdxList=zeros(numVertices,1);
totalInSetIdxList=0;
numSets=0;

for v=1:numVertices
    if(~reached(v))
        depthFirstSearch(v);
    end
end

%Shrink to the number of sets.
setStartList=setStartList(1:numSets);
%SetIdxList will include all of the vertices and thus does not need to be
%resized.

%Record the lengths of the sets to be returned.
setLengths=zeros(numSets,1);
if(numSets==1)
    setLengths=numVertices;
else
    for curSet=1:(numSets-1)
        setLengths(curSet)=setStartList(curSet+1)-setStartList(curSet);
    end
    setLengths(numSets)=numVertices-setStartList(numSets)+1;
end

if(nargout>3)%If desired, actually return the subgraphs.
    subgraphList=cell(numSets,1);
    
    for curSubgraph=1:numSets
        numSubgraphNodes=setLengths(curSubgraph);
        adjSubMat=zeros(numSubgraphNodes,numSubgraphNodes);
        
        %Fill in the values of the adjacency submatrix.
        baseIdx=setStartList(curSubgraph);
        subGraphIdx=setIdxList(baseIdx:(baseIdx+numSubgraphNodes-1));
        for v=1:numSubgraphNodes
            for w=1:numSubgraphNodes
                adjSubMat(v,w)=adjMat(subGraphIdx(v),subGraphIdx(w));
            end
        end
        subgraphList{curSubgraph}=adjSubMat;
    end
end


function depthFirstSearch(v)
    count1=count1+1;
    dfsNum(v)=count1;
    reached(v)=true;

    %Push v onto unfinished
    numInUnfinished=numInUnfinished+1;
    unfinished(numInUnfinished)=v;
    inUnfinished(v)=true;

    %Push v onto roots
    numInRoots=numInRoots+1;
    roots(numInRoots)=v;

    %Go through all of the edges involving v
    for w=1:numVertices
        %If no edge exists.
        if(adjMat(v,w)==Inf)
            continue;
        end

        if(~reached(w))
            depthFirstSearch(w);
        elseif(inUnfinished(w))%Now, we merge components
            while(dfsNum(roots(numInRoots))>dfsNum(w))
                %pop roots
                numInRoots=numInRoots-1;
            end
        end
    end
    
    if(v==roots(numInRoots))
        %This is added to keep track of where the current strongly
        %connected subgraph begins.
        totalInSetIdxList=totalInSetIdxList+1;
        setIdxList(totalInSetIdxList)=v;
        numSets=numSets+1;
        setStartList(numSets)=totalInSetIdxList;
        
        while(1)
            w=unfinished(numInUnfinished);
            numInUnfinished=numInUnfinished-1;
            inUnfinished(w)=false;
            %w is an element of the strongly connected set with root v.
            if(v==w)
                break;
            else
                totalInSetIdxList=totalInSetIdxList+1;
                setIdxList(totalInSetIdxList)=w;
            end
        end
        %pop roots
        numInRoots=numInRoots-1;
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
