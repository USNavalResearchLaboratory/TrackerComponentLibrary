function [xi,w]=secondOrderCubPointsBingham(A,lambda)
%%SECONDORDERCUBPOINTSBINGHAM Generate second-order cubature points for
%                   numeric integrals involving the Bingham distribution,
%                   which is defined on a unit (hyper-)sphere. The Bingham
%                   distribution is a central version of the Fisher-Bingham
%                   distribution. The PDF is
%                   f(x|A,g)=(1/c)*exp(-x'*A*x+g'*x)
%                   with g=0 and c a normalizing constant.
%
%INPUTS: A A numDimXnumDim type of inverse covariance matrix for the
%          distribution. It must be symmetric. It does not have to be
%          positive definite.
%   lambda An optional parameter that affects how much weights is assigned
%          to a point placed at the pole. This must be between 0 and 1.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%The algorithm implementes is essentially that given in Section III of [1],
%but corrected. An assumption in the algorithm of the paper is that the 
%
%REFERENCES:
%[1] I. Gilitschenski, G. Kurz, S. J. Julier, and U. D. Hanebeck,
%    "Unscented orientation estimation based on the Bingham distribution,"
%    IEEE Transactions on Automatic Control, accepted 2015.
%   
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(lambda))
    lambda=0.5;
end

numDim=size(A,1);

%First, find the omega terms. These are expressed in terms of the
%eigenvalues of the distribution. The definition in Section II of [1]
%requires that the eigenvalues be sorted in ascending order. That does not
%matter here.
[M,Z]=eig(A);

omega=diag(FisherBinghamD.covApprox(Z));

xi = zeros(numDim,2*numDim-1);
w = zeros(2*numDim-1,1);

%The sample at the mean.
xi(numDim,1)=1;
for k=1:(numDim-1)
    p=omega(k)+(1-lambda)*(omega(numDim)/(numDim-1));
    alpha=asin(sqrt(omega(k)/p));
    xi(end,2*k)=cos(alpha);
    xi(end,2*k+1)=cos(alpha);
    xi(k,2*k)=sin(alpha);
    xi(k,2*k+1)=-sin(alpha);

    w(2*k)=p/2;
    w(2*k+1)=p/2;
end
w(1)=1-sum(w(2:end));

%The points generated are mirrored to get all of the points.
xi=[xi,-xi];
w=[w;w]/2;

w=w/sum(w);%Normalize.

xi=M*xi;

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
