function [adjMat, distMat] = genRipsGraph(data, thresh)
%%GENRIPSGRAPH Given a set of N data points in k-dimensional space, a Rips
%              graph is generated using the given threshold parameter. A
%              Rips graph is an adjacency matrix, where the diagonal of the
%              matrix is zeros (something is not adjacent with itself).
%
%INPUTS:  data A k x N matrix with each column representing a data point
%              in k-dimensional space.
%       thresh A user-defined threshold that determines which data points
%              are considered "neighbors". neighbors have an l2 norm <=
%              this distance.
%
%OUTPUTS:  adjMat An NXN adjacency matrix for the resulting Rips graph.
%         distMat An NXN  matrix whose (i,j)th coordinate is the pairwise
%                 distance between nodes i and j.
%
%EXAMPLE 1:
%This is just a simple example with a set of 1D data indicating that the
%adjacency matrix produced by this function agrees with the actual
%adjacency matrix that one can form by inspecting the data.
% thresh=1;
% data=[0,10,2,3,8,7,100,101,1];
% adjMat=genRipsGraph(data, thresh);
% adjMatTrue=zeros(9,9);
% adjMatTrue(3,4)=1;
% adjMatTrue(4,3)=1;
% adjMatTrue(5,6)=1;
% adjMatTrue(6,5)=1;
% adjMatTrue(7,8)=1;
% adjMatTrue(8,7)=1;
% adjMatTrue([1,3],9)=1;
% adjMatTrue(9,[1,3])=1;
% allEqual=all(adjMatTrue(:)==adjMat(:))
%
%EXAMPLE 2:
%This example demonstrates the results of applying this function to a pair
%of nonconvex clusters. Nonzero entries are displayed for several
%thresholds to demonstrate the growth in connections corresponding to the
%threshold parameter.
% numPoints = 2e4;
% points = FibonacciGrid(numPoints,3); % last input is numdim+1
% points(1:2,:) = 5*points(1:2,:)-2.5; % rescale to be roughly centered on origin
% mu_meas1 = [0;.5];
% dst_meas1 = 1;
% std_meas1 = .2;
%
% % Creating two nonconvex densities to sample
% mu_prior1 = mu_meas1;
% mu_prior1(1) = -dst_meas1;
% C_prior1 = [1^2 0; 0 1^2];
%
% prior1 = @(x,y) exp(-1/2 * dot(([x;y]-mu_prior1),(C_prior1\([x;y]-mu_prior1)),1) );
% meas1  = @(x,y) (sum(([x;y]-mu_meas1).^2,1));
% ll1    = @(x,y) exp(-1/2 * ((meas1(x,y)-dst_meas1)/std_meas1).^2 );
% post1  = @(x,y) prior1(x,y) .* ll1(x,y);
%
% prior2 = @(x,y) exp(-1/2 * dot(([x;y]+mu_prior1),(C_prior1\([x;y]+mu_prior1)),1) );
% meas2  = @(x,y) (sum(([x;y]+mu_meas1).^2,1));
% ll2    = @(x,y) exp(-1/2 * ((meas2(x,y)-dst_meas1)/std_meas1).^2 );
% post2  = @(x,y) prior2(x,y) .* ll2(x,y);
%
% f = post1(points(1,:),points(2,:)) + post2(points(1,:),points(2,:));
%
% keep = points(3,:)<f; % sample points below the density
%
% % Plot the points and their likelihood
% figure(1)
% scatter3(points(1,keep),points(2,keep),f(keep),'ok','LineWidth',2,'CData',f(keep))
% title({sprintf("Rejection Sampling of Banana Density using %0.e Fibonacci Points",numPoints),
%       sprintf("%d Points Were Kept",sum(keep))})
% figure(2)
% tr = delaunay(points(1,keep),points(2,keep));
% trimesh(tr,points(1,keep),points(2,keep),f(keep))
%
% % Save the kept points as in a variable
% % sampled points in a numdim+1-by-numsamples matrix
% % row 1 is x coordinate, row 2 is y coordinate, row 3 is likelihood value
% points = points(:,keep);
% fkeep = f(keep);
%
% % Generating the Rips Graph
% for thresh = [0.2,0.3,0.4,0.5,0.6]
%    [adjMat, distMat] = genRipsGraph(points,thresh);
%    figure
%    spy(adjMat)
% end
%
%July 2024 Michael Kardos, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Stores number of data points to avoid repeated computation
dataLen = size(data,2);

%Initializes pairwise distance matrix
distMat = zeros(dataLen,dataLen);

for i = 1:dataLen
    for j = i+1:dataLen
        %Stores Euclidean distance between data points i and j to the
        %(i,j)th entry of distMat
        distMat(i,j) = norm(data(:,i) - data(:,j));
        distMat(j,i) = distMat(i,j);
    end
end

%Joins vertices whose pairwise distances is below the set threshold
adjMat = distMat <= thresh;
adjMat = adjMat - eye(dataLen);
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
