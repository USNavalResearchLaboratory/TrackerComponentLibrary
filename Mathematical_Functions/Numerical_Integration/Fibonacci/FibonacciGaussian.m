function [points,weights] = FibonacciGaussian(N,mu,P)
%%FIBONACCIGAUSSIAN Creates a grid of points based on a transformation of a
%                   Fibonacci grid as given in [1].
%
%INPUT:
% N: The number of points to sample.
% mu: The D-by-1 vector indicating the Gaussian mean.
% P: The D-by-D covariance matrix. This must be positive definite.
%
%OUTPUT:
% points: A D-by-N matrix containing sample points for the given Gaussian.
% weights: A collection of weights representing the Gaussian PDF at each
%          sample point.
%
%EXAMPLE 1: Generates a collection of sample points and their respective
%           Gaussian weights in 1D and plots them.
% mu = 5;
% P = 2;
% [points, weights] = FibonacciGaussian(100,mu,P);
% stem(points(1,:),weights,'k',"filled","MarkerEdgeColor",'b')
%
%EXAMPLE 2: Generates a collection of sample points and their respective
%         Gaussian weights in 2D and plots them.
% mu = [5;3];
% P = [2,-0.7;-0.7,1];
% [points, weights] = FibonacciGaussian(100,mu,P);
% scatter(points(1,:),points(2,:),'k',"filled","MarkerEdgeColor",'b',"CData",weights)
% hold on
% drawEllipse(mu,P\eye(2),[],'r',"LineWidth",2)
% hold off
%
%EXAMPLE 3: Generates a collection of sample points and their respective
%           Gaussian weights in 3D and plots them.
% mu = [5;3;4];
% P = [2,-0.7,0.1;-0.7,1,-0.3;0.1,-0.3,2];
% [points, weights] = FibonacciGaussian(100,mu,P);
% scatter3(points(1,:),points(2,:),points(3,:),'k',"filled","MarkerEdgeColor",'b',"CData",1e3*weights)
% hold on
% drawEllipse(mu,P\eye(3),[],"FaceAlpha",0.1,"EdgeAlpha",0.1)
% hold off
%
%REFERENCE:
%[1] D. Frisch and U. D. Hanebeck, "Deterministic gaussian sampling with
%    generalized fibonacci grids," in 2021 IEEE 24th International
%    Conference on Information Fusion (FUSION), IEEE, 2021, pp. 1-8.
%
%December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
D = size(mu,1);
points = FibonacciGrid(N,D);
points = sqrt(2)*erfinv(2*points-1);
nu = sqrt(sum(points.^2,2)/N);
for idx = 1:D
    points(idx,:) = points(idx,:)/nu(idx);
end
[V,E] = eig(P);
points = V*sqrt(abs(E))*points+mu;
weights = GaussianD.PDF(points,mu,P);
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