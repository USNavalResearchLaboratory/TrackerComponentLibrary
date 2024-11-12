function lattice = FibonacciSphereLattice(N,alg,algFG,center)
%%FIBONACCISPHERELATTICE Generates 3D vectors corresponding to the
%                        spherical Fibonacci lattice.
%
%INPUT:
% N: An integer which determines how many points to generate in the
%    lattice.
% alg: An integer indicating which Fibonacci point set to generate:
%      0) (Default) A direct generation method (doesn't call to
%         FibonacciGrid) which uses a golden ratio spiral based on [1] and
%         [2]. This method only requires the input parameter N and is not
%         affected by the other inputs.
%      1) A Lambert cylindrical equal-area transform from the unit square
%         to the sphere. This generates the unit square Fibonacci points
%         using FibonacciGrid. The input parameters algFG and center 
%         determine which generation method FibonacciGrid uses. This is
%         described in [3].
%      2) A spherical coordinates (angle only) to Cartesian transform from
%         the unit square expanded to [0,pi)-by-[0,2*pi). This uses the
%         spher2Cart function for the transformation. This generates the
%         unit square Fibonacci points using FibonacciGrid. The input
%         parameters algFG and center determine which generation method
%         FibonacciGrid uses.
% algFG: An integer indicating which Fibonacci point set to generate. See
%        FibonacciGrid for more information. If not given or empty,
%        algFG=0. Note that algFG=1 returns the next Fibonacci number of
%        points greater than N. The default is 0.
% center: A Boolean which is true if the set should be centered and false
%         otherwise. This only affects options 1 and 2 of algFG. The
%         default is false.
%
%OUTPUT:
% lattice: A 3-by-N matrix containing column vectors corresponding to each
%          point in the lattice.
%
%EXAMPLE: Generates a lattice of 100 points and plots the points and convex
%         hull.
% lat = FibonacciSphereLattice(100);
% figure(1); clf
% plot3(lat(1,:),lat(2,:),lat(3,:),'o')
% grid on
% figure(2); clf
% trisurf(convhull(lat(1,:),lat(2,:),lat(3,:)),lat(1,:),lat(2,:),lat(3,:))
%
%REFERENCES:
%[1] B. Keinert, M. Innmann, M. Sainger, and M. Stamminger, "Spherical
%    Fibonacci mapping," ACM Trans. Graph., vol. 34, no. 6, 2015,
%    issn: 0730-0301. doi: 10.1145/2816795.2818131. [Online].
%    Available: https://doi.org/10.1145/2816795.2818131.
%[2] A. Gonzalez, "Measurement of areas on a sphere using Fibonacci and
%    latitude-longitude lattices," Mathematical Geosciences, vol. 42,
%    no. 1, pp. 49-64, 2010.
%[3] C Aistleitner, J. S. Brauchart, and J Dick, "Point sets on the sphere
%    S2 with small spherical cap discrepancy," Discrete & Computational
%    Geometry, vol. 48, pp. 990-1024, 4 2012, issn: 1432-0444.
%    doi: 10.1007/s00454-012-9451-3. [Online].
%    Available: https://doi.org/10.1007/s00454-012-9451-3.
%
%March 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
if ~exist('alg','var') || isempty(alg)
    alg = 0;
end
if ~exist('algFG','var') || isempty(algFG)
    algFG = 0;
end
if ~exist('center','var') || isempty(center)
    center = false;
end
switch alg
    case 0 % A direct generation method from [1] and [2]
        lattice = zeros(3,N);
        gratio = (1+sqrt(5))/2;
        for idx = 0:N-1
            theta = 2*pi*idx/gratio;
            phi = acos(1-2*(idx+0.5)/N);
            x = sin(phi)*cos(theta);
            y = sin(theta)*sin(phi);
            z = cos(phi);
            lattice(:,idx+1) = [x;y;z];
        end
    case 1 % A Lambert cylindrical equal-area mapping from [3]
        ugrid = FibonacciGrid(N,2,algFG,center);
        px = ugrid(1,:);
        py = ugrid(2,:);
        lattice = [2*sqrt(px-px.^2).*cos(2*pi*py);
            2*sqrt(px-px.^2).*sin(2*pi*py);
            1-2*px];
    case 2 % A spherical coordinates to Cartesian coordinates conversion
        ugrid = FibonacciGrid(N,2,algFG,center);
        px = pi*ugrid(1,:);
        py = 2*pi*ugrid(2,:);
        lattice = spher2Cart([px;py]);
    otherwise
        error("Unknown algorithm choice.")
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