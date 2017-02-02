function [xi,w]=secondOrderSpherCubPoints(numDim)
%%SECONDORDERSPHERCUBPOINTS Generate second-order cubature points for
%               integration over a unit spherical (or hyperspherical)
%               region (weighting function is just 1).
%
%INPUTS: numDim  An integer specifying the dimensionality of the points
%                to be generated.
%
%OUTPUTS:   xi      A numDim X numCubaturePoints matrix containing the
%                   cubature points. (Each "point" is a vector)
%           w       A numCubaturePoints X 1 vector of the weights
%                   associated with the cubature points.
%
%This implements formula Sn 2-1 in [1], pg. 267, numDim+1points.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

V=(2/numDim)*pi^(numDim/2)/gamma(numDim/2);

w=V/(numDim+1)*ones(numDim+1,1);

xi=zeros(numDim,numDim+1);

i=0:numDim;
for k=1:fix(numDim/2)
    xi(2*k-1,i+1)=sqrt(2/(numDim+2))*cos(2*i*k*pi/(numDim+1));
    xi(2*k,i+1)=sqrt(2/(numDim+2))*sin(2*i*k*pi/(numDim+1));
end

if(mod(numDim,2)==1)
    xi(numDim,:)=(-1).^i/sqrt(numDim+2);
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
