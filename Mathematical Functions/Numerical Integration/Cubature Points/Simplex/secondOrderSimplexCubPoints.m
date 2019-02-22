function [xi,w]=secondOrderSimplexCubPoints(numDim)
%%SECONDORDERSIMPLEXCUBPOINTS Generate second-order cubature points for
%               integration over the simplex of points containing all
%               points such that sum(x)<=1 and all x>=0, where x is a
%               numDimX1 vector.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points. (Each "point" is a vector)
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%This implements formula T_n 2-2 in [1], pg. 307, (numDim+1)*(numDim+2)/2
%points
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

n=numDim;
V=1/factorial(n);
r=1/2;
B=V*(2-n)/((n+1)*(n+2));
C=V*4/((n+1)*(n+2));

numPoints=(numDim+1)*(numDim+2)/2;
xi=zeros(numDim,numPoints);
w=zeros(numPoints,1);

xiExtended=genAllMultisetPermutations([zeros(n,1);1]);
num2Add=size(xiExtended,2);
xi(:,1:num2Add)=xiExtended(1:n,:);
w(1:num2Add)=B;
curStart=num2Add+1;

xiExtended=genAllMultisetPermutations([zeros(n-1,1);r;r]);
xi(:,curStart:end)=xiExtended(1:n,:);
w(curStart:end)=C;

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
