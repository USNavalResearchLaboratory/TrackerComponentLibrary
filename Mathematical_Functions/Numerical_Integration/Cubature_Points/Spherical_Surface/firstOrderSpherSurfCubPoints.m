function [xi,w]=firstOrderSpherSurfCubPoints(numDim)
%%FIRSTORDERSPHERSURFCUBPOINTS Obtain third-order cubature points for
%                   integration over the surface of the unit sphere
%                   (weighting function is just 1).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>=1.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%The formula for the points is the simple formula on page 431 of [1].
%
%REFERENCES:
%[1] A. Genz and J. Monahan, "Stochastic integration rules for infinite
%    regions," SIAM Journal on Scientific Computing, vol. 19, no. 2, pp.
%    426-439, Mar. 1998.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

absUm=2*pi^(numDim/2)/gamma(numDim/2);

w=repmat(absUm/2,[2,1]);
xi=zeros(numDim,2);
xi(1,1)=1;
xi(1,2)=-1;

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
