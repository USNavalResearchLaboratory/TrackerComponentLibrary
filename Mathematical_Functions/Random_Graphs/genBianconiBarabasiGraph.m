function [A,availableNodes] = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode,fitness,initA,omitSelfLoops,omitMultiLoops,delProb,isDirected,maxInDegree,maxOutDegree)
%%GENBIANCONIBARABASIGRAPH Generates a random instance of the
%                          Bianconi-Barabasi model using the linearized
%                          chord diagram assumptions given in [1]. Note the
%                          Barabasi-Albert model of [2] can be generated
%                          as a special case of this function.
%
%INPUT:
% finalN: Number of nodes in the final graph.
% initN: Number of nodes in the initialized graph. If initA is omitted, the
%        initial nodes each begin with a single self-loop.
% addedEdgesPerNode: Number of edges to be added from each additional node.
%                    This must be a nonnegative number. Note this is the
%                    number of proposed edges. Other settings such as
%                    omitMultiLoops or maxInDegree may result in fewer
%                    edges connected to a node. Defaults to finalN if
%                    omitted, empty, or Inf.
% fitness: An optional length finalN vector indicating relative fitness
%          between the nodes in terms of favorability for receiving edges
%          from new nodes. If omitted or an empty vector is given, this
%          is a vector of ones, resulting in the Barabasi-Albert model.
% initA: An optional initN-by-initN weighted adjacency matrix which is used
%        to initialize edges between the first initN nodes. If omitted, the
%        implementation is equivalent to passing 2*eye(initN) if
%        isDirected is false and eye(initN) if isDirected is true. A test
%        for symmetry is performed to ensure only bidirectional edges are
%        used.
% omitSelfLoops: A boolean which, if true, will cause the generated graph
%                to ignore self-loops. The default if omitted or empty is
%                false, consistent with the original Bianconi-Barabasi
%                model. Note this does not modify initA.
% omitMultiLoops: A boolean which if true, will prevent multiple loops the
%                 undirected case (this is impossible in the directed
%                 case). The default if omitted or empty is false,
%                 consistent with the original Bianconi-Barabasi model.
% delProb: A number between [0,1] determining the frequency with which a
%          random node is deleted from the network on average. The default
%          if omitted or empty is to never delete nodes.
% isDirected: A boolean indicating whether the graph to be generated should
%             be directed or undirected. Undirected only checks for
%             connections between nodes once. Directed checks for both the
%             incoming connection and the outgoing connection on all nodes.
%             The default if omitted or an empty array is passed is false.
% maxInDegree: A non-negative integer specifying the maximum number of
%              allowed incoming edges for each node. For an undirected
%              graph, this is used as the degree and maxOutDegree is
%              ignored. The default if omitted or an empty array is passed
%              is Inf.
% maxOutDegree: A non-negative integer specifying the maximum number of
%               allowed incoming edges for each node. For an undirected
%               graph, this parameter is ignored. The default if omitted or
%               an empty array is passed is Inf.
%
%OUTPUT:
% A: The weighted adjacency matrix for the final graph with deleted nodes
%    included. If finalN<1e3, this is just a matrix. If finalN>=1e3, a
%    sparse matrix is returned. To retrieve the adjacency matrix for just
%    the remaining (not deleted) nodes, use
%    A(availableNodes,availableNodes).
% availableNodes: A vector of node indices 1:finalN including all nodes
%                 which were not deleted.
%
%NOTE: The definition of degree used here counts a bidirectional edge
%      twice, once for the outgoing node and once for the incomming node.
%      Therefore, self loops are represented by a weight of 2 and edges
%      between two distinct nodes split that weight between the (i,j) and
%      (j,i) elements of A, each gaining a weight of 1.
%
%EXAMPLE 1: Generate a Barabasi-Albert model using the default fitness and
%           initialized graph.
% finalN = 30;
% initN = 3;
% addedEdgesPerNode = 2;
% A = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode);
% g = graph(A);
% plot(g,"NodeCData",degree(g))
% colorbar
% title("Degree of Nodes in Bianconi-Barabasi Model")
%
%EXAMPLE 2: Generate a Bianconi-Barabasi model using a specified fitness
%           vector which only allows connections to odd nodes.
% finalN = 30;
% initN = 3;
% addedEdgesPerNode = 2;
% fitness = zeros(finalN,1);
% fitness(1:2:finalN) = 1;
% A = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode,fitness);
% g = graph(A);
% plot(g,"NodeCData",degree(g))
% colorbar
% title("Degree of Nodes in Bianconi-Barabasi Model")
%
%EXAMPLE 3: Generate a Bianconi-Barabasi model using a specified fitness
%           vector which only allows connections to odd nodes and a
%           complete graph on the set of initialized nodes.
% finalN = 30;
% initN = 3;
% addedEdgesPerNode = 2;
% fitness = zeros(finalN,1);
% fitness(1:2:finalN) = 1;
% Ain = eye(initN);
% Ain = Ain+ones(size(Ain));
% gin = graph(Ain);
% figure(1)
% plot(gin)
% A = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode,fitness,Ain);
% g = graph(A);
% figure(2)
% plot(g,"NodeCData",degree(g))
% colorbar
% title("Degree of Nodes in Bianconi-Barabasi Model")
%
%EXAMPLE 4: Similar to example 3, but this time including a nonzero
%           probability of deleting nodes. Note that modifying the number
%           of nodes may require modifying the multiplier on the plotting
%           marker size in order to render the graph presentable. A graph
%           displaying the degree distribution is also produced.
% finalN = 1e3;
% initN = 10;
% addedEdgesPerNode = 2;
% fitness = zeros(finalN,1);
% fitness(1:2:finalN) = 1;
% Ain = eye(initN);
% Ain = Ain+ones(size(Ain));
% gin = graph(Ain);
% figure(1)
% plot(gin)
% title("Initial Graph")
% [A,nodeIdxs] = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode,fitness,Ain,false,true,0.01);
% g = graph(A);
% figure(2)
% markerscaling = 2e2;
% plot(g,"NodeCData",degree(g),"Layout","force")
% colorbar
% title("Degree of Nodes in Bianconi-Barabasi Model")
% figure(3)
% [GC,GV] = groupcounts(degree(g));
% semilogy(GV,GC,'bo',"MarkerFaceColor",'b')
% xlabel("Degree")
% ylabel("# of Nodes (log scale)")
% title("Number of Nodes with a Given Degree")
%
%EXAMPLE 5: An instance of the restricted, directed Bianconi-Barabasi model
%           is generated using an initial graph of 3 nodes which have
%           self-loops and 2 edges each. The remaining nodes are added
%           without self-loops and with a maximum of 1 outgoing edge. A
%           figure displaying the initial graph is produces along with
%           figures showing the indegree and outdegree for each node. Note
%           that this function removes self-loops, but does not force the
%           initial graph to meet the outdegree restriction.
% finalN = 30;
% initN = 3;
% addedEdgesPerNode = 2;
% fitness = zeros(finalN,1);
% fitness(1:2:finalN) = 1;
% Ain = eye(initN);
% Ain = Ain+ones(size(Ain));
% gin = digraph(Ain);
% figure(1)
% plot(gin)
% title("Initial Graph")
% A = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode,fitness,Ain,true,true,[],true,Inf,1);
% g = digraph(A);
% figure(2)
% plot(g,"NodeCData",indegree(g))
% colorbar
% title("Indegree of Nodes in Restricted, Directed Bianconi-Barabasi Model")
%
%REFERENCES:
%[1] G. Bianconi and A.-L Barabasi, "Competition and multiscaling in
%    evolving networks," Europhysics Letters (EPL), vol. 54, no. 4,
%    pp. 436-442, 2001.
%[2] A.-L. Barabasi and R. Albert, "Emergence of scaling in random
%    networks," Science, vol. 286, no. 5439, pp. 509-512, 1999.
%
%August 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if finalN<initN
    error('initN must be less than final N')
end
if nargin<3 || isempty(addedEdgesPerNode) || isinf(addedEdgesPerNode)
    addedEdgesPerNode = finalN;
end
% If fitness is not given, reduce to a Barabasi-Albert model.
if nargin<4 || isempty(fitness)
    fitness = ones(finalN,1);
else
    try
        fitness = reshape(fitness,finalN,1);
    catch
        error('fitness vector must have N elements.')
    end
end
if nargin<9 || isempty(isDirected)
    isDirected = false;
end
if nargin<5 || isempty(initA)
    if isDirected
        initA = eye(initN);
    else
        initA = 2*eye(initN);
    end
else
    if ~issymmetric(initA)
        error('initA must be a symmetric matrix.')
    elseif initN~=size(initA,1)
        error('initA must be of shape initN-by-initN.')
    end
end
if nargin<6 || isempty(omitSelfLoops)
    omitSelfLoops = false;
end
if nargin<7 || isempty(omitMultiLoops)
    omitMultiLoops = false;
end
if nargin<8 || isempty(delProb)
    delProb = 0;
elseif delProb<0 || delProb>1
    error('delRate must be on the interval [0,1]')
end
if nargin<10 || isempty(maxInDegree)
    maxInDegree = Inf;
elseif maxInDegree<0
    error("maxInDegree must be nonnegative.")
end
if nargin<11 || isempty(maxOutDegree)
    maxOutDegree = Inf;
elseif maxOutDegree<0
    error("maxOutDegree must be nonnegative.")
end
if ~isDirected
    maxOutDegree = maxInDegree;
end

if finalN<1e3
    A = zeros(finalN);
else
    A = sparse(finalN,finalN);
end
A(1:initN,1:initN) = initA;
availableNodes = true(1,initN);
for i = initN+1:finalN
    % First delete a node with rate delRate
    delNode = rand()<delProb;
    if delNode
        idx = find(availableNodes,randi(sum(availableNodes)));
        idx = idx(end);
        node2Del = availableNodes(idx);
        A(node2Del,:) = 0;
        A(:,node2Del) = 0;
        availableNodes(node2Del) = 0;
        fitness(node2Del) = 0;
    end
    
    % Then add new node and edges.
    r = rand(1,addedEdgesPerNode);
    % Preference is a product of node degree or indegree and fitness
    preference = sum(A(1:i-1,1:i-1),1)'.*fitness(1:i-1);
    [plink,I] = sort([preference;fitness(i)]/(fitness(i)+sum(preference)));
    plink = cumsum(plink);
    for j = 1:length(r)
        link = find(r(j)<plink,1);
        if ~omitMultiLoops || (omitMultiLoops && A(i,I(link))==0)
            if I(link)~=i
                if sum(A(i,:),2)<maxOutDegree && sum(A(:,I(link)),1)<maxInDegree
                    if ~isDirected
                        A(I(link),i) = A(I(link),i)+1;
                    end
                    A(i,I(link)) = A(i,I(link))+1;
                end
            else
                if ~omitSelfLoops && sum(A(i,:),2)<maxOutDegree && sum(A(:,I(link)),1)<maxInDegree
                    if ~isDirected
                        A(I(link),i) = A(I(link),i)+2;
                    else
                        A(I(link),i) = A(I(link),i)+1;
                    end
                end
            end
        end
    end
    availableNodes(i) = 1;
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
