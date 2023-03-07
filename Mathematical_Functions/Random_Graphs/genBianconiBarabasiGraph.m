function [A,availableNodes] = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode,fitness,initA,delProb)
%%GENBIANCONIBARABASIGRAPH Generates a random instance of the
%                          Bianconi-Barabási model using the linearized
%                          chord diagram assumptions given in [1]. Note the
%                          Barabási-Albert model of [2] can be generated
%                          as a special case of this function.
%
%INPUT:
% finalN: Number of nodes in the final graph.
% initN: Number of nodes in the initialized graph. If initA is omitted, the
%        initial nodes each begin with a single self-loop.
% addedEdgesPerNode: Number of edges to be added from each additional node.
% fitness: An optional length finalN vector indicating relative fitness
%          between the nodes in terms of favorability for receiving edges
%          from new nodes. If omitted or an empty vector is given, this
%          is a vector of ones, resulting in the Barabási-Albert model.
% initA: An optional initN-by-initN weighted adjacency matrix which is used
%        to initialize edges between the first initN nodes. If omitted, the
%        implementation is equivalent to passing 2*eye(initN). A test for
%        symmetry is performed to ensure only bidirectional edges are used.
% delProb: A number between [0,1] determining the frequency with which a
%          random node is deleted from the network on average.
%
%OUTPUT:
% A: The weighted adjacency matrix for the final graph with deleted nodes
%    omitted. If finalN<1e3, this is just a matrix. If finalN>=1e3, a
%    sparse matrix is returned.
% availableNodes: A vector of node indices 1:finalN including all nodes
%                 which were not deleted.
%
%NOTE: The definition of degree used here counts a bidirectional edge
%      twice, once for the outgoing node and once for the incomming node.
%      Therefore, self loops are represented by a weight of 2 and edges
%      between two distinct nodes split that weight between the (i,j) and
%      (j,i) elements of A, each gaining a weight of 1.
%
%EXAMPLE 1: Generate a Barabási-Albert model using the default fitness and
%           initialized graph.
% finalN = 30;
% initN = 3;
% addedEdgesPerNode = 2;
% A = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode);
% g = graph(A);
% plot(g,"MarkerSize",degree(g))
%
%EXAMPLE 2: Generate a Bianconi-Barabási model using a specified fitness
%           vector which only allows connections to odd nodes.
% finalN = 30;
% initN = 3;
% addedEdgesPerNode = 2;
% fitness = zeros(finalN,1);
% fitness(1:2:finalN) = 1;
% A = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode,fitness);
% g = graph(A);
% plot(g,"MarkerSize",degree(g))
%
%EXAMPLE 3: Generate a Bianconi-Barabási model using a specified fitness
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
% plot(gin)
% A = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode,fitness,Ain);
% g = graph(A);
% figure
% plot(g,"MarkerSize",degree(g))
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
% plot(gin)
% [A,nodeIdxs] = genBianconiBarabasiGraph(finalN,initN,addedEdgesPerNode,fitness,Ain,0.01);
% g = graph(A);
% figure
% markerscaling = 2e2;
% plot(g,"MarkerSize",markerscaling*(degree(g)+1)/sum(degree(g)),"Layout","force")
% figure
% [GC,GV] = groupcounts(degree(g));
% semilogy(GV,GC,'bo',"MarkerFaceColor",'b')
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

% If initial graph adjacency matrix is not provided, set to default.
% Otherwise, check input for symmetry and dimensions.
if nargin<5
    initA = 2*eye(initN);
else
    if ~issymmetric(initA)
        error('initA must be a symmetric matrix.')
    elseif initN~=size(initA,1)
        error('initA must be of shape initN-by-initN.')
    end
end

% If deletion rate is not provided, do not delete nodes.
% Otherwise, check that the provided rate is in [0,1].
if nargin<6
    delProb = 0;
elseif delProb<0 || delProb>1
    error('delRate must be on the interval [0,1]')
end

if finalN<1e3
    A = zeros(finalN);
else
    A = sparse(finalN,finalN);
end
% Initialize first nodes
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
    % Preference is a product of node degree and fitness
    preference = sum(A(1:i-1,1:i-1),2).*fitness(1:i-1);
    [plink,I] = sort([preference;fitness(i)]/(fitness(i)+sum(preference)));
    plink = cumsum(plink);
    for j = 1:length(r)
        link = find(r(j)<plink,1);
        if I(link)~=i
            A(I(link),i) = A(I(link),i)+1;
            A(i,I(link)) = A(i,I(link))+1;
        else
            A(I(link),i) = A(I(link),i)+2;
        end
    end
    availableNodes(i) = 1;
end
A = A(availableNodes,availableNodes);
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
