function [x,z]=constrainedLSEq(A,b,C,d)
%%CONSTRAINEDLSEQ Find x to minimize norm(A*x-b) under the constraint that
%                 C*x=d. It is assumed that the (m+p)Xn stacked matrix
%                 [A;C] has linearly independent columns (left invertible)
%                 and C has linearly independent rows (right invertible).
%                 This means that p<=n<=m+p. This is equality constrained
%                 least squares.
%
%INPUTS: A A real mXn matrix
%        b A real mX1 vector.
%        C A real pXn matrix.
%        d A real pX1 vector.
%
%OUTPUTS: x The solution to the equality constrained least squares problem.
%           An mX1 vector.
%         z The Lagrange multipliers used in the optimization.
%
%In Chapter 14 of [1], it is shown that x solves the constrained least
%squares problem if and only if there exists a z such that
%[A'*A, C'; C, zeros(p,p)]*[x;z]=[A'*b;d]
%This is a set of n+p equations and unknowns. Thus, this function just
%solves the above system of equations. z is a set of Lagrange multipliers.
%
%REFERENCES:
%[1] S. Boyd and L. Vandenberghe, Vectors, Matrices, and Least Squares,
%    Sep. 2015, draft version of book. [Online].
%    Available: http://www.seas.ucla.edu/~vandenbe/133A/133A-textbook.pdf
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(A,2);

p=size(C,1);
H=A'*A;
c=A'*b;
xz=[H,C';C,zeros(p,p)]\[c;d];

x=xz(1:m);
z=xz((m+1):end);

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
