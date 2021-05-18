function v=polyIntValNewton(z,a,c,K)
%%POLYINTVALNEWTON Evaluate integrals of a set of polynomial given in
%                  Newton's form at a desired set of points.
%
%INPUTS: z A 1XnumZ or 2XnumZ matrix of the scalar points at which the
%          integrals should be evaluated. If z is a 2XnumZ matrix, then
%          this function will return the definite integrals (for each of
%          the numDim equations) between each pair of points (from z(1,:)
%          to z(2,:).
%        a A numDimXn matrix of polynomial coefficients for numDim
%          polynomials, as described below.
%        c A numDim X(n-1) matrix of the control points associated with
%          each of the numDim polynomial in Newton's form. When performing
%          Hermite interpolation using a polynomial returned by the
%          HermiteInterpPoly function, c can be a matrix holding the
%          values where the interpolating polynomial matches the data for
%          numDim calls to the HermiteInterpPoly function (one for each
%          dimension).
%        K The optional scalar constant of integration to be added to the
%          evaluated values. Assumes K=0 if not supplied.
%
%OUTPUTS: pInt A numDimXnumZ matrix of the integrals for each of the numZ
%              points or pairs of points for eahc of the numDim equations.
%
%A scalar polynomial function in Newton's form evaluated at point z has the
%form
%y(z)=a(1)+sum_{k=1}^{n-1}a(k+1)(z-c(1))*(z-c(2))*...*(z-c(k))
%This function will evaluate numDim integrals of such functions
%
%The underlying approach is based on VALUE in Chapter 19 of [1], which
%provides a recursion for evaluating polynomials in Newton form. The
%integral algorithm comes from integrating the recursion. The evaluation of
%the integral has the form
%Y(z)=sum{k=1}^{n}a(k)sum_{m=1}^{k}(-1)^(k-m)*z^m/m*sum(prod(C_{k-m}(c(1:k-1))
%where C_{k}(c) is all the possible combinations of the vector c, taken k
%at a time, and prod(C) is the product of all the elements in each set.
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%
%February 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(K))
    K=0;
end

numZ=size(z,2);
numDim=size(a,1);

v=zeros(numDim,numZ);
for curDim=1:numDim
    v(curDim,:)=polyIntValNewtonScalar(z,a(curDim,:),c(curDim,:),K);
end

end

function v=polyIntValNewtonScalar(z,a,c,K)

n=length(a);

v=ones(size(z))*K;
for k=1:n
    kterm=0;
    for l=1:k
        C=nchoosek(c(1:k-1),k-l);
        if isempty(C)
            C=1;
        end
        lterm=(-1)^(k-l)*z.^l/l*sum(prod(C,2));
        kterm=kterm+lterm;
    end
    v=v+a(k)*kterm;
end

%If definite integrals are desired.
if size(z,1)==2
    v=v(2,:)-v(1,:);
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
