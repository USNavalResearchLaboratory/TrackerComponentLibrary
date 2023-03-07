function A = genErdosRenyiGraph(N,pConnect)
%%GENERDOSRENYIGRAPH Generates an instance of the Erdös-Renyi model
%                    described in [1-3].
%
%INPUT:
% N: The number of nodes used for generating the graph.
% pConnect: The probability of any chosen edge being realized.
%
%OUTPUT:
% A: The adjacency matrix for the final graph.
%
%EXAMPLE: Generates an instance of the Erdös-Renyi model.
% N = 30;
% pConnect = 0.3;
% A = genErdosRenyiGraph(N,rewireProb);
% g = graph(A);
% plot(g,"MarkerSize",degree(g))
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

A = zeros(N);
for i = 1:N-1
    for j = i+1:N
        A(i,j) = rand()<pConnect;
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
