function A = genErdosRenyiGraph(N,pConnect,isDirected,maxInDegree,maxOutDegree)
%%GENERDOSRENYIGRAPH Generates an instance of the Erdos-Renyi model
%                    described in [1-3].
%
%INPUT:
% N: The number of nodes used for generating the graph.
% pConnect: The probability of any chosen edge being realized.
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
% A: The adjacency matrix for the final graph.
%
%EXAMPLE 1: Generates an instance of the conventional Erdos-Renyi model.
% N = 30;
% pConnect = 0.3;
% A = genErdosRenyiGraph(N,pConnect);
% g = graph(A);
% plot(g,"NodeCData",degree(g))
% colorbar
% title("Degree of Nodes in Erdos-Renyi Model")
%
%EXAMPLE 2: Generates an instance of the Erdos-Renyi digraph model.
% N = 30;
% pConnect = 0.1;
% A = genErdosRenyiGraph(N,pConnect,true);
% g = digraph(A);
% plot(g,"NodeCData",indegree(g))
% colorbar
% title("Indegree of Nodes in Directed Erdos-Renyi Model")
%
%EXAMPLE 3: Generates an instance of the restricted, directed Eros-Renyi
%           model.
% N = 30;
% pConnect = 0.3;
% A = genErdosRenyiGraph(N,pConnect,true,5,2);
% g = digraph(A);
% plot(g,"NodeCData",indegree(g))
% colorbar
% title("Indegree of Nodes in Restricted, Directed Erdos-Renyi Model")
%
%REFERENCES:
%[1] P. Erdos, A. Renyi, et al., "On the evolution of random graphs,"
%    Publ. Math. Inst. Hung. Acad. Sci, vol. 5, no. 1, pp. 17-60, 1960.
%[2] P. Erdos and A. Renyi, "On random graphs i," Publ. Math. Debrecen,
%    vol. 6, pp. 290-297, 1959.
%[3] E. N. Gilbert, "Random graphs," The Annals of Mathematical
%    Statistics, vol. 30, no. 4, pp. 1141-1144, 1959.
%
%August 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if pConnect>1 || pConnect<0
    error('pConnect must be in the interval [0,1]')
end
if nargin<3 || isempty(isDirected)
    isDirected = false;
end
if nargin<4 || isempty(maxInDegree)
    maxInDegree = Inf;
elseif maxInDegree<0
    error("maxInDegree must be nonnegative.")
end
if nargin<5 || isempty(maxOutDegree)
    maxOutDegree = Inf;
elseif maxOutDegree<0
    error("maxOutDegree must be nonnegative.")
end
if ~isDirected
    maxOutDegree = maxInDegree;
end

A = zeros(N);
if ~isDirected
    rows = 1:N-1;
    rows = rows(randperm(length(rows)));
    for i = rows
        cols = i+1:N;
        cols = cols(randperm(length(cols)));
        for j = cols
            if sum(A(i,:),2)<maxOutDegree && sum(A(:,j),1)<maxInDegree
                A(i,j) = rand()<pConnect;
            end
        end
    end
    A = A+A';
else
    rows = 1:N;
    rows = rows(randperm(length(rows)));
    for i = rows
        cols = setdiff(1:N,i);
        cols = cols(randperm(length(cols)));
        for j = cols
            if sum(A(i,:),2)<maxOutDegree && sum(A(:,j),1)<maxInDegree
                A(i,j) = rand()<pConnect;
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
