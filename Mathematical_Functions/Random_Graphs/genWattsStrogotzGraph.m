function A = genWattsStrogotzGraph(N,K,rewireProb)
%%GENWATTSSTROGOTZGRAPH Generates a random instance of the Watts-Strogotz
%                       model as described in [1].
%
%INPUT:
% N: The number of nodes to be used in generating the graph.
% K: The number of connected neighbors for each node. This must be an even
%    number.
% rewireProb: The probability that any single edge connecting a given node
%             to its neighbor will be changed to connect to a different
%             node chosen uniformly at random from the set of nodes which
%             do not at that time share an edge with the given node.
%
%OUTPUT:
% A: The adjacency matrix for the final graph.
%
%EXAMPLE: Generates an instance of the Watts-Strogotz model.
% N = 30;
% K = 4;
% rewireProb = 0.3;
% A = genWattsStrogotzGraph(N,K,rewireProb);
% g = graph(A);
% plot(g,"MarkerSize",degree(g),"Layout","circle")
%
%REFERENCES:
%[1] D. J. Watts and S. H. Strogatz, "Collective dynamics of 'small-world'
%    networks," Nature, vol. 393, no. 6684, pp. 440-442, 1998.
%
%August 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if mod(K,2)==1
    error('K must be even')
end
if rewireProb>1 || rewireProb<0
    error('rewireProb must be in the interval [0,1]')
end

A = zeros(N);
for i = 1:N-1
    for j = i+1:N
        A(i,j) = mod(abs(i-1-(j-1)),N-K/2)<=K/2;
    end
end
for i = 1:N
    for j = i+1:i+K/2
        jNode = mod(j-1,N)+1;
        if rand()<rewireProb
            A(i,jNode) = 0;
            free = find(A(i,:)==0);
            newNode = i;
            while newNode==i
                newNode = free(randi(length(free)));
            end
            if newNode<i
                A(newNode,i) = 1;
            else
                A(i,newNode) = 1;
            end
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