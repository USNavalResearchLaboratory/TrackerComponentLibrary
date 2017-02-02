function [xi,w]=spherSurfPoints2GaussPoints(xiSurf,wSurf,order,beta)
%%SPHERSURFPOINTS2GAUSSPOINTS Given cubature points of a certain order for
%               integration over the surface of a unit (hyper-)sphere,
%               obtain cubature points for integration over a multivariate
%               Gaussian PDF times |x|^beta (beta=0 means just the Gaussian
%               PDF). THe weighting function is thus
%               w(x)=1/(2*pi)^(numDim/2)*norm(x)^beta*exp(-x'*x/2).
%
%INPUTS: xiSurf, wSurf A numDimXnumPoints set of cubture points and the
%               corresponding numPointsX1 set of cubature weights for
%               integration over the surface of the unit sphere.
%         order The positive order of the spherical surface cubature points
%               passed. This will also be the order of spherical points
%               produced by the algorithm.
%          beta An integer specifying the exponent of the norm(x) term in
%               the weighting function (the thing times the normal PDF in
%               the weighting function). beta>-numDim. If omitted or an
%               empty matrix is passed, beta=0 is used.
%
%The method of taking a spherical shell and a formula for a 1D integral
%over exp(-x^2) times |x|^beta and transforming it for multiple dimensions
%is described in Chapter 2.8 of [1]. The transformation to a Gaussian PDF
%is performed by scaling the points and weights.
%
%As an example, consider obtaining fourteenth-order points for integrating
%over a Gaussian PDF in 3D:
% [xiSurf,wSurf]=fourteenthOrderSpherSurfCubPoints(3);
% [xi,w]=spherSurfPoints2GaussPoints(xiSurf,wSurf,14);
% %As a test for the points, consider the integral of z^14:
% sum(bsxfun(@times,sum(xi(1,:).^14,1),w'))
%One will find that the result is the same as the 14th moment in the first
%dimension given by GaussianD.momentGenFun(zeros(3,1),eye(3),[14;0;0])
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(beta))
   beta=0; 
end

numDim=size(xiSurf,1);
deltaStep=size(xiSurf,2);

num1DPoints=ceil((order+1)/2);
[r,A]=quadraturePoints1D(num1DPoints,12,numDim-1+beta);
%We discard the negative half of the points as mentioned in Chapter 2.8.
%That changes the points from integrating from -1 to 1 to integrating from
%0 to 1. This only works, because of the symmetry of the points about 0.
sel=(r>=0);
r=r(sel);
A=A(sel);
num1DPoints=length(r);

%Allocate space
numPoints=num1DPoints*deltaStep;
xi=zeros(numDim,numPoints);
w=zeros(numPoints,1);

curStart=1;
for i=1:num1DPoints
    xi(:,curStart:(curStart+deltaStep-1))=r(i)*xiSurf;
    w(curStart:(curStart+deltaStep-1))=A(i)*wSurf;
    curStart=curStart+deltaStep;
end

%We now have cubature points and weights for integrating with a weighting
%function of w(x)=|x|^beta*exp(-x'*x). We will transform these to
%integrating over the normal o-I distribution times |x|^beta.
w=pi^(-numDim/2)*sqrt(2)^(beta)*w;
xi=sqrt(2)*xi;
    
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
