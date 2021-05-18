function [xi,w]=secondOrderCrossPolyCubPoints(numDim)
%%SECONDORDERCROSSPOLYCUBPOINTS Generate second order cubature points for
%               integration over the region sum(abs(x))<=1 with weighting
%               function w(x)=1. In two dimensions, this is a square, in
%               three an octohedron, in n a cross polytope (or
%               n-dimensional octahedron).
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%Algorithm G_n 2-1 in [1], pg. 303, with n+1 points is used.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

V=exp(numDim*log(2)-gammaln(numDim+1));
n=numDim;

w=V/(1+n)*ones(n+1,1);

xi=zeros(n,n+1);
i=0:n;
for k=1:fix(n/2)
    xi(2*k-1,i+1)=sqrt(4/((n+1)*(n+2)))*cos(2*i*k*pi/(n+1));
    xi(2*k,i+1)=sqrt(4/((n+1)*(n+2)))*sin(2*i*k*pi/(n+1));
end

if(mod(n,2)==1)
    xi(end,i+1)=(-1).^i.*sqrt(2/((n+2)*(n+2)));
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
