function Xfib = FibonacciGrid(N,D)
%%FIBONACCIGRID Generates a Fibonacci grid on the unit volume in D
%               dimensions using N points.
%
%INPUT:
% N: The number of points to sample.
% D: The dimension of the space to be sampled.
%
%OUTPUT:
% Xfib: A D-by-N matrix containing points of a Fibonacci grid.
%
%EXAMPLE: A scatter plot of the 3D Fibonacci grid.
% points = FibonacciGrid(100,3);
% scatter3(points(1,:),points(2,:),points(3,:),'filled')
%
%REFERENCE:
%[1] D. Frisch and U. D. Hanebeck, "Deterministic gaussian sampling with
%    generalized fibonacci grids," in 2021 IEEE 24th International
%    Conference on Information Fusion (FUSION), IEEE, 2021, pp. 1-8.
%[2] R. J. Purser, "Generalized fibonacci grids; a new class of structured,
%    smoothly adaptive multi-dimensional computational lattices," Office
%    note (National Centers for Environmental Prediction (U.S.)) 2008.
%    [Online]. Available: https://repository.library.noaa.gov/view/noaa/6956.
%
%December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
Vd = FibEigen(D);
s = sum(abs(Vd),2);
shc = max(s);
N0 = ceil(N^(1/D));
delta = 1/N0;
N1 = ceil(shc/delta)+2;
if mod(N,2)~=mod(N1,2)
    N1 = N1+1;
end
r = (1:N1)'*delta; 
r = r - sum(r/N1,"all");
Xregcell = cell(1,D);
[Xregcell{:}] = ndgrid(r);
for idx = 1:D
    Xregcell{idx} = Xregcell{idx}(:)';
end
Xreg = cat(1,Xregcell{:});
Xrot = Vd'*Xreg;
valid = all(Xrot>-0.5 & Xrot<0.5,1);
Xfib = Xrot(:,valid);
N2 = size(Xfib,2);
Xfib = sortrows(Xfib',1)';
u = (N2-N)/2;
Xfib = Xfib(:,u+1:N2-u);
b = 0.5-0.5/N;
m = max(Xfib,[],2);
for idx = 1:D
    Xfib(idx,:) = (Xfib(idx,:)/m(idx))*b;
end
Xfib = Xfib+ones(size(Xfib))/2;
end

function Vd = FibEigen(d)
%%FIBEIGEN Generates a normalized eigenvector matrix for the Fibonacci
%          grid as given in [1] and [2].
%
%INPUT:
% d: Dimension in which the lattice will be generated.
%
%OUTPUT:
% Vd: A d-by-d matrix of normalized columns of eigenvectors.
%
%REFERENCES:
%[1] D. Frisch and U. D. Hanebeck, "Deterministic gaussian sampling with
%    generalized fibonacci grids," in 2021 IEEE 24th International
%    Conference on Information Fusion (FUSION), IEEE, 2021, pp. 1-8.
%[2] R. J. Purser, "Generalized fibonacci grids; a new class of structured,
%    smoothly adaptive multi-dimensional computational lattices," Office
%    note (National Centers for Environmental Prediction (U.S.)) 2008.
%    [Online]. Available: https://repository.library.noaa.gov/view/noaa/6956.
%
%December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
Vd = zeros(d);
for i = 1:d
    for j = 1:d
        Vd(i,j) = cos((2*i-1)*(2*j-1)*pi/(4*d+2));
    end
end
Vd = normc(Vd);
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