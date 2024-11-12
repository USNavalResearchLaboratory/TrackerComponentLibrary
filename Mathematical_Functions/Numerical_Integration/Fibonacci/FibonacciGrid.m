function Xfib = FibonacciGrid(N,D,alg,center)
%%FIBONACCIGRID Generates a Fibonacci grid on the unit volume in D
%               dimensions using N points.
%
%INPUT:
% N: The number of points to sample.
% D: The dimension of the space to be sampled.
% alg: An integer indicating which Fibonacci point set to generate:
%      0) (Default) A grid technique due to Purser and based on the
%         algorithm given in [1].
%      1) The standard Fibonacci points which have been proven to have the
%         lowest discrepancy possible in dimension D=2. Noted in [2] and
%         [3]. The number of points output will be the nearest Fibonacci
%         number greater than N. This algorithm is only defined for D=2.
%      2) A variant based on the algorithm in option 1, but loses
%         periodicity in one dimension for non-Fibonacci number of points
%         N. Also noted in [3]. This algorithm is only defined for D=2.
% center: A Boolean which is true if the set should be centered and false
%         otherwise. The standard point sets for algorithms 1 and 2 have a
%         point at (0,0) for every N. If this parameter is true, then the
%         lattice is shifted by 1/(2*N) in both dimensions to provide a
%         lattice centered within the unit square. This only affects
%         options 1 and 2 of alg. The default is false.
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
%[3] D. Frisch and U. D. Hanebeck, "Deterministic von mises–fisher
%    sampling on the sphere using fibonacci lattices," in 2023 IEEE
%    Symposium Sensor Data Fusion and International Conference on
%    Multisensor Fusion and Integration (SDF-MFI), IEEE, Nov. 2023,
%    pp. 1-8, isbn: 979-8-3503-8258-7. doi: 10.1109/SDF-MFI59545.2023.10361396.
%    [Online]. Available: https://ieeexplore.ieee.org/document/10361396/.
%
%December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
if ~exist('alg','var') || isempty(alg)
    alg = 0;
end
if ~exist('center','var') || isempty(center)
    center = false;
end
if alg>0 && D~=2
    error("Algorithms 1 and 2 are only defined for 2 dimensions.")
end

switch alg
    case 0 % Swinbank-Purser Lattice
        Vd = FibEigen(D);
        %Get bounds on box containing rotated points
        s = sum(abs(Vd),2);
        shc = max(s);
        
        %Compute the uniform sample vector (1d)
        N0 = ceil(N^(1/D));
        N1 = ceil(shc*N0); %[1] has a +2 here, but not needed
        if mod(N,2)~=mod(N1,2)
            N1 = N1+1;
        end
        r = (1:N1)'/N0; 
        r = r - sum(r/N1,"all");
        
        %Get a D-dimensional grid sampled at r in each dimension
        Xregcell = cell(1,D);
        [Xregcell{:}] = ndgrid(r);
        for idx = 1:D
            Xregcell{idx} = Xregcell{idx}(:)';
        end
        Xreg = cat(1,Xregcell{:});
        
        %Rotate grid and prune points in d=2,3,...,D
        Xrot = Vd'*Xreg;
        valid = all(Xrot>-0.5 & Xrot<0.5,1);
        Xfib = Xrot(:,valid);
        
        %Decide how many other points need to be pruned, and prune them in d=1
        N2 = size(Xfib,2);
        Xfib = sortrows(Xfib',1)';
        u = (N2-N)/2; % number of remaining points to be removed (one side)
        
        %Check if too few points. If so, redo point generation.
        if u<0
            %Compute the uniform sample vector (1d)
            N0 = ceil(N^(1/D)+0.5);
            N1 = ceil(shc*N0); %[1] has a +2 here, but not needed
            if mod(N,2)~=mod(N1,2)
                N1 = N1+1;
            end
            r = (1:N1)'/N0; 
            r = r - sum(r/N1,"all");
            
            %Get a D-dimensional grid sampled at r in each dimension
            Xregcell = cell(1,D);
            [Xregcell{:}] = ndgrid(r);
            for idx = 1:D
                Xregcell{idx} = Xregcell{idx}(:)';
            end
            Xreg = cat(1,Xregcell{:});
            
            %Rotate grid and prune points in d=2,3,...,D
            Xrot = Vd'*Xreg;
            valid = all(Xrot>-0.5 & Xrot<0.5,1);
            Xfib = Xrot(:,valid);
            
            %Decide how many other points need to be pruned, and prune them in d=1
            N2 = size(Xfib,2);
            Xfib = sortrows(Xfib',1)';
            u = (N2-N)/2; % number of remaining points to be removed
        end
        Xfib = Xfib(:,u+1:N2-u);
        
        %Rescale to unit volume
        b = 0.5-0.5/N;
        m = max(Xfib,[],2);
        for idx = 1:D
            Xfib(idx,:) = (Xfib(idx,:)/m(idx))*b;
        end
        Xfib = Xfib+ones(size(Xfib))/2;
    case 1 % Niederreiter-Sloan Lattice
        k = FibonacciNumInv(N,2); % get the index of the next highest Fib num
        Fk = FibonacciNum(k);
        Fk1 = FibonacciNum(k-1);
        Xfib = zeros(D,Fk);
        fibvec = [1/Fk;Fk1/Fk];
        for idx = 0:Fk-1
            Xfib(:,idx+1) = (idx*fibvec);
        end
        if center
            Xfib = Xfib+1/(2*Fk);
        end
        Xfib = mod(Xfib,1);
    case 2 % Kronecker-Jacob Lattice
        Xfib = zeros(D,N);
        fibvec = [1/N;(sqrt(5)-1)/2];
        for idx = 0:N-1
            Xfib(:,idx+1) = (idx*fibvec);
        end
        if center
            Xfib = Xfib+1/(2*N);
        end
        Xfib = mod(Xfib,1);
    otherwise
        error("Unknown algorithm choice.")
end
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