function [xi,w]=fourthOrderSimplexCubPoints(numDim)
%%FOURTHORDERSIMPLEXCUBPOINTS Generate fourth-order cubature points for
%               integration over the simplex of points containing all
%               points such that sum(x)<=1 and all x>=0, where x is a
%               numDimX1 vector.
%
%INPUTS: numDim An integer specifying the dimensionality of the points to
%               be generated. numDim>=3.
%
%OUTPUTS: xi A numDim X numCubaturePoints matrix containing the cubature
%            points (Each "point" is a vector).
%          w A numCubaturePoints X 1 vector of the weights associated with
%            the cubature points.
%
%Formula T_n 4-1 in [1], pg. 311,
%factorial(numDim+4)/(factorial(4)*factorial(numDim) points, is used.
%
%REFERENCES:
%[1] A.H. Stroud, Approximate Calculation of Multiple Integrals. Cliffs,
%    NJ: Prentice-Hall, Inc., 1971.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

if(numDim<3)
   error('This algorithm requires that numDim>=3') 
end

n=numDim;
V=1/factorial(n);

r=1/4;
s=3/4;
t=1/2;
B1=V*(-3*n^3+17*n^2-58*n+72)/(3*(n+1)*(n+2)*(n+3)*(n+4));
B2=V*16*(n^2-5*n+12)/(3*(n+1)*(n+2)*(n+3)*(n+4));
B3=V*4*(n^2-9*n+12)/((n+1)*(n+2)*(n+3)*(n+4));
B4=V*64*(4-n)/(2*(n+1)*(n+2)*(n+3)*(n+4));
B5=V*256/((n+1)*(n+2)*(n+3)*(n+4));

numPoints=fix(exp(gammaln(n+4+1)-(gammaln(4+1)+gammaln(n+1))));
xi=zeros(n,numPoints);
w=zeros(numPoints,1);
curStart=1;

xiExtended=genAllMultisetPermutations([zeros(n,1);1]);
num2Add=size(xiExtended,2);
xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
w(curStart:(curStart+num2Add-1))=B1;
curStart=curStart+num2Add;

xiExtended=genAllMultisetPermutations([zeros(n-1,1);r;s]);
num2Add=size(xiExtended,2);
xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
w(curStart:(curStart+num2Add-1))=B2;
curStart=curStart+num2Add;

xiExtended=genAllMultisetPermutations([zeros(n-1,1);t;t]);
num2Add=size(xiExtended,2);
xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
w(curStart:(curStart+num2Add-1))=B3;
curStart=curStart+num2Add;

xiExtended=genAllMultisetPermutations([zeros(n-2,1);r;r;t]);
num2Add=size(xiExtended,2);
xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
w(curStart:(curStart+num2Add-1))=B4;
curStart=curStart+num2Add;

xiExtended=genAllMultisetPermutations([zeros(n-3,1);r;r;r;r]);
num2Add=size(xiExtended,2);
xi(:,curStart:(curStart+num2Add-1))=xiExtended(1:n,:);
w(curStart:(curStart+num2Add-1))=B5;

%Get rid of zero-weight points, which will occur if n==4.
sel=~(w==0);
w=w(sel);
xi=xi(:,sel);

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
