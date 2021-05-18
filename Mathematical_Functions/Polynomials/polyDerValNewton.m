function [pDer,p]=polyDerValNewton(z,a,c,numDer)
%%POLYDERVALNEWTON Evaluate derivatives of a set of polynomials given in
%                  Newton's form at a desired set of points.
%
%INPUTS: z A numZX1 or 1XnumZ vector of the points at which the derivatives
%          should be evaluated.
%        a A numDimXn matrix of polynomial coefficients for numDim
%          polynomials, as described below.
%        c A numDimX(n-1) matrix of the control points associated with
%          each of the numDim polynomial in Newton's form. When performing
%          Hermite interpolation using a polynomial returned by the
%          HermiteInterpPoly function, c can be a matrix holding the
%          values where the interpolating polynomial matches the data for
%          numDim calls to the HermiteInterpPoly function (one for each
%          dimension).
%   numDer The positive integer number of derivatives to return. If this
%          parameter is omitted, then numDer=1 is used.
%
%OUTPUTS: pDer A numDimXnumDerXnumZ hypermatrix of the derivatives for each
%              of the numZ points for each of the numDim dimensions.
%            p A numDimXnumZ vector of the function values at the points.
%
%A polynomial function in Newton's form evaluated at point z has the form
%y(z)=a(1)+sum_{k=1}^{n-1}a(k+1)(z-c(1))*(z-c(2))*...*(z-c(k))
%
%The algorithm is based on VALUE in Chapter 19 of [1], which provides a
%recursion for evaluating polynomials in Newton form. The derivative
%algorithm comes from differentiating the recursion.
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    numDer=1;
end

numZ=length(z);
numDim=size(a,1);

%Allocate space
pDer=zeros(numDim,numDer,numZ);
p=zeros(numDim,numZ);

for curDim=1:numDim
    [pDer(curDim,:,:),p(curDim,:)]=polyDerValNewtonScalar(z,a(curDim,:),c(curDim,:),numDer);
end
end

function [pDer,p]=polyDerValNewtonScalar(z,a,c,numDer)

n=length(a);
numZ=length(z);

%Make sure that z is a column vector.
z=z(:)';

pDerList=zeros(numDer+1,numZ);
pDerList(1,:)=a(n);
for k=(n-1):-1:1
    for curDir=(numDer+1):-1:2
        pDerList(curDir,:)=(curDir-1)*pDerList(curDir-1,:)+(z-c(k)).*pDerList(curDir,:);
    end
    
    pDerList(1,:)=(z-c(k)).*pDerList(1,:)+a(k);
end

p=pDerList(1,:);
pDer=pDerList(2:end,:);
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
