function [xi,w]=fourteenthOrderCubPoints(numDim,beta,randomize)
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
%     randomize If this parameter is true, then the points will be
%               multiplied by a random orthonormal rotation matrix. This
%               does not change the moments up to the order of the points.
%               This randomization is done in [1] and [2] to lessen various
%               effects that arise when using points in the same
%               orientation repeatedly in tracking. The default if this
%               parameter is omitted or an empty matrix is passed is false.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%This function gets fourteenth-order spherical surface points using
%fourteenthOrderSpherSurfCubPoints
%and then transforms them to work with a multidimensional Gaussian PDF
%using
%spherSurfPoints2GaussPoints.
%
%REFERENCES:
%[1] O. Straka, D. Duník, and M. Simandl, "Randomized unscented Kalman
%    filter in tracking," in Proceedings of the 15th International
%    Conference on Information Fusion, Singapore, 9-12 Jul. 2012, pp.
%    503-510.
%[2] J. Duník, O. Straka, and M. Simandl, "The development of a randomised
%    unscented Kalman filter," in Proceedings of the 18th World Congress,
%    The International Federation of Automatic Control, Milan, Italy, 28
%    Aug. - 2 Sep. 2011, pp. 8-13.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(numDim~=3)
   error('Only 3D points are supported'); 
end

if(nargin<3||isempty(randomize))
    randomize=false;
end

if(nargin<2||isempty(beta))
    beta=0;
end

[xiSurf,wSurf]=fourteenthOrderSpherSurfCubPoints(3);
[xi,w]=spherSurfPoints2GaussPoints(xiSurf,wSurf,14,beta);

if(randomize)
    R=randOrthoMat(numDim);

    numPoints=length(w);
    for curPoint=1:numPoints
        xi(:,curPoint)=R*xi(:,curPoint);
    end
end
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
