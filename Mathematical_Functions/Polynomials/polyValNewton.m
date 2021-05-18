function v=polyValNewton(z,a,c)
%POLYVALNEWTON Evaluate a set of polynomial given in Newton's form at a
%              desired set of points. This is useful when given a set of
%              Hermite interpolating polynomials for interpolating a
%              vector.
%
%INPUTS: z A numPointsX1 or 1XnumPoints vector of scalar points where the
%          values of the numDim polynomials are desired.
%        a A numDimXn matrix of polynomial coefficients for numDim
%          polynomials, as described below.
%        c A numDim X(n-1) matrix of the control points associated with
%          each of the numDim polynomial in Newton's form, which might be
%          the control points used in the HermiteInterpPoly function.
%
%OUTPUTS: v The numDimXnumPoints matrix of values of the polynomials
%           evaluated at the points in z.
%
%A polynomial function in Newton's form evaluated at point z has the form
%y(z)=a(1)+sum_{k=1}^{n-1}a(k+1)(z-c(1))*(z-c(2))*...*(z-c(k))
%This function just evaluates multiple scalar polynomials.
%
%The algorithm is based on  VALUE from Chapter 19 of [1] to handle vector
%values; a loop is used here.
%
%Newton's polynomial form arises when dealing with Hermite interpolating
%polynomials as is discussed in Section 2.1.3 of [2] and in Chapters 3.3
%and 3.4 of [3].
%
%REFERENCES:
%[1] A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%    and Calculators, 2nd ed. New York: Academic press, 1978.
%[2] J. Stoer and R. Bulirsch, Introduction to Numerical Analysis, 3rd ed.
%    New York: Springer, 2002.
%[3] R. L. Burden and J. D. Faires, Numerical Analysis, 9th ed. Boston, MA:
%    Brooks/Cole, 2011.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(a,1);
numPoints=length(z);

v=zeros(numDim,numPoints);
for curDim=1:numDim
   v(curDim,:)=polyValNewtonScalar(z,a(curDim,:),c(curDim,:));
end
end

function v=polyValNewtonScalar(z,a,c)
%Below is VALUE from Chapter 19 of [1]
%
%REFERENCES:
%A. Nijenhuis and H. S. Wilf, Combinatorial Algorithms for Computers
%and Calculators, 2nd ed. New York: Academic press, 1978.

n=length(a);

v=ones(size(z))*a(n);
for k=(n-1):-1:1
    v=(z-c(k)).*v+a(k);
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
