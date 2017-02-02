function xi=transformCubPoints(xi,mu,S)
%%TRANSFORMCUBPOINTS Transform a set of cubature points for a normal(0,I)
%                    distribution to approximation a normal distribution
%                    with mean mu and covariance matrix S*S'.
%
%INPUTS:    xi      A numDim X numCubaturePoints matrix containing the
%                   cubature points for a normal(0,I) distribution. (Each
%                   "point" is a vector). The weights associated with the
%                   points do not need to be transformed when applying the
%                   points to other Gaussian distributions
%           mu      A numDim X 1 vector that is the desired mean of the
%                   transformed points.
%           S       A numDum X numDim lower-triangular square root of the
%                   covariance matrix of the desired distribution. If P is 
%                   the covariance matrix of the desired distribution, then
%                   one can obtain S as S=chol(P,'lower');
%
%OUTPUTS:   xi      The cubature points transformed to be useful in
%                   approximation integrals over a normal distribution with 
%                   mean mu and covariance matrix S*S'.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xi=bsxfun(@plus,S*xi,mu);
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
