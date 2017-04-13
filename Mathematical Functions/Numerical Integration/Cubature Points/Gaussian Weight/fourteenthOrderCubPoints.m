function [xi,w]=fourteenthOrderCubPoints(numDim,beta)
%%FOURTEENTHORDERCUBPOINTS Generate fourteenth order cubature points for
%           integration  over a multivariate Gaussian PDF times |x|^beta
%           (beta=0 means just the Gaussian PDF). The weighting function is
%           thus w(x)=1/(2*pi)^(numDim/2)*norm(x)^beta*exp(-x'*x/2).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. Currently, only 3D points are supported.
%          beta An integer specifying the exponent of the norm(x) term in
%               the weighting function (the thing times the normal PDF in
%               the weighting function). beta>-numDim. If omitted or an
%               empty matrix is passed, beta=0 is used.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%This function gets fourteenth-order spherical surface points using
%fourteenthOrderSpherSurfCubPoints
%and then transforms them to work with a multidimensional gausian PDF using
%spherSurfPoints2GaussPoints.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(numDim~=3)
   error('Only 3D points are supported'); 
end

if(nargin<2||isempty(beta))
    beta=0;
end

[xiSurf,wSurf]=fourteenthOrderSpherSurfCubPoints(3);
[xi,w]=spherSurfPoints2GaussPoints(xiSurf,wSurf,14,beta);
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
