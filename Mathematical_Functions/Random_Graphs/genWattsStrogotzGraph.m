function A = genWattsStrogotzGraph(N,K,rewireProb,isDirected,maxInDegree)
%%GENWATTSSTROGOTZGRAPH Generates a random instance of the Watts-Strogotz
%                       model as described in [1].
%
%INPUT:
% N: The number of nodes to be used in generating the graph.
% K: The number of connected neighbors for each node. This must be an even
%    number. If isDirected is true, this is the outdegree of each node. If
%    K>N, then K is set to the largest even number less than or equal to N.
% rewireProb: The probability that any single edge connecting a given node
%             to its neighbor will be changed to connect to a different
%             node chosen uniformly at random from the set of nodes which
%             do not at that time share an edge with the given node.
% isDirected: A boolean indicating whether the graph to be generated should
%             be directed or undirected. Undirected only checks for
%             connections between nodes once. Directed checks for both the
%             incoming connection and the outgoing connection on all nodes.
%             The default if omitted or an empty array is passed is false.
% maxInDegree: A non-negative integer specifying the maximum number of
%              allowed incoming edges for each node. The default if omitted
%              or an empty array is passed is Inf.
%
%OUTPUT:
% A: The adjacency matrix for the final graph.
%
%EXAMPLE 1: Generates an instance of the Watts-Strogotz model.
% N = 30;
% K = 4;
% rewireProb = 0.3;
% A = genWattsStrogotzGraph(N,K,rewireProb);
% g = graph(A);
% plot(g,"NodeCData",degree(g),"Layout","circle")
% colorbar
% title("Degree of Nodes in Watts-Strogatz Model")
%
%EXAMPLE 2: Generates an instance of the directed Watts-Strogotz model.
% N = 30;
% K = 2;
% rewireProb = 0.3;
% A = genWattsStrogotzGraph(N,K,rewireProb,true);
% g = digraph(A);
% plot(g,"NodeCData",indegree(g),"Layout","circle")
% colorbar
% title("Indegree of Nodes in Directed Watts-Strogatz Model")
%
%EXAMPLE 3: Generates an instance of the restricted, directed
%           Watts-Strogotz model.
% N = 30;
% K = 4;
% rewireProb = 0.3;
% A = genWattsStrogotzGraph(N,K,rewireProb,true,5);
% g = digraph(A);
% plot(g,"NodeCData",indegree(g),"Layout","circle")
% colorbar
% title("Indegree of Nodes in Restricted, Directed Watts-Strogatz Model")
%
%REFERENCES:
%[1] D. J. Watts and S. H. Strogatz, "Collective dynamics of 'small-world'
%    networks," Nature, vol. 393, no. 6684, pp. 440-442, 1998.
%
%August 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if mod(K,2)==1 || K<0
    error('K must be a nonnegative even number.')
end
if K>N
    if mod(N,2)==0
        K = N;
    else
        K = N-1;
    end
end
if rewireProb>1 || rewireProb<0
    error('rewireProb must be in the interval [0,1]')
end
if nargin<4 || isempty(isDirected)
    isDirected = false;
end
if nargin<5 || isempty(maxInDegree)
    maxInDegree = Inf;
elseif maxInDegree<0
    error("maxInDegree must be nonnegative.")
end

A = zeros(N);
if ~isDirected
    rows = 1:N;
    rows = rows(randperm(length(rows)));
    for i = rows
        cols = setdiff(1:N,i);
        cols = cols(randperm(length(cols)));
        cols = cols(mod(abs(i-1-(cols-1)),N-K/2)<=K/2);
        for j = cols
            if sum(A(:,j),1)<maxInDegree && sum(A(:,i),1)<maxInDegree
                A(i,j) = 1;
                A(j,i) = 1;
            end
        end
    end

    rows = 1:N;
    rows = rows(randperm(length(rows)));
    for i = rows
        cols = setdiff(i+1:i+K/2,i);
        cols = cols(randperm(length(cols)));
        for j = cols
            jNode = mod(j-1,N)+1;
            if rand()<rewireProb
                A(i,jNode) = 0;
                A(jNode,i) = 0;
                free = find(A(i,:)==0);
                newNode = i;
                while newNode==i
                    newNode = free(randi(length(free)));
                end
                if sum(A(:,newNode),1)<maxInDegree
                    A(i,newNode) = 1;
                    A(newNode,i) = 1;
                else
                    A(i,jNode) = 1;
                    A(jNode,i) = 1;
                end
            end
        end
    end
else
    rows = 1:N;
    rows = rows(randperm(length(rows)));
    for i = rows
        cols = setdiff(1:N,i);
        cols = cols(randperm(length(cols)));
        cols = cols(mod(abs(i-1-(cols-1)),N-K/2)<=K/2);
        for j = cols
            if sum(A(:,j),1)<maxInDegree
                A(i,j) = 1;
            end
        end
    end

    rows = 1:N;
    rows = rows(randperm(length(rows)));
    for i = rows
        cols = setdiff(i-K/2:i+K/2,i);
        cols = cols(randperm(length(cols)));
        for j = cols
            jNode = mod(j-1,N)+1;
            if rand()<rewireProb
                A(i,jNode) = 0;
                free = find(A(i,:)==0);
                newNode = i;
                while newNode==i
                    newNode = free(randi(length(free)));
                end
                if sum(A(:,newNode),1)<maxInDegree
                    A(i,newNode) = 1;
                else
                    A(i,jNode) = 1;
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