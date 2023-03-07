function A = probBasedRandomGraph(N,P)
%PROBBASEDRANDOMGRAPH Generates a graph with random connections whose
%                     probability is based on pair-wise probabilities.
%
%INPUTS:
% N: The number of nodes in the graph.
% P: A N-by-N matrix containing probabilities of pair-wise connections. If
%    an undirected graph is desired, this should be an upper triangular
%    matrix. Otherwise, a directed graph is assumed.
%
%OUTPUT:
% A: The N-by-N adjacency matrix for the realized graph. This has 2 on the
%    diagonal for self-loops to count both ends of the edges.
%
%EXAMPLE 1: Generate a random digraph using randomly generated pairwise
%           probabilities.
% rng(42)
% P = rand(10);
% A = probBasedRandomGraph(10,P);
% g = digraph(A);
% plot(g,"MarkerSize",indegree(g))
%
%EXAMPLE 2: Generate a random undirected graph using randomly generated
%           pairwise probabilities.
% rng(42)
% P = triu(rand(10));
% A = probBasedRandomGraph(10,P);
% g = graph(A);
% plot(g,"MarkerSize",degree(g))
%
%August 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if any(P>1 | P<0)
    error('Edge probabilities must be in [0,1]')
end

if istriu(P)
    directed = false;
elseif istril(P)
    P = P';
    directed = false;
else
    directed = true;
end

r = rand(N);
A = zeros(N);
if directed
    for i = 1:N
        for j = 1:N
            if r(i,j)<P(i,j)
                if i~=j
                    A(i,j) = 1;
                else
                    A(i,j) = 2;
                end
            end
        end
    end
else %undirected
    for i = 1:N
        for j = i:N
            if r(i,j)<P(i,j)
                A(i,j) = 1;
            end
        end
    end
    A = A+A';
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
