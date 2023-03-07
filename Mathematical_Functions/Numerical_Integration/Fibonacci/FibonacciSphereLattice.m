function lattice = FibonacciSphereLattice(N)
%%FIBONACCISPHERELATTICE Generates 3D vectors corresponding to the
%                        spherical Fibonacci lattice.
%
%INPUT:
% N: An integer which determines how many points to generate in the
%    lattice.
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
%
%March 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

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