function [x,M]=semidefQuadProgEqConst(B,a,A,b)
%%SEMIDEFQUADPROGEQCONST Solve a quadratic programming problem involving a
%       semidefinite matrix (or a positive defintie matrix) where there are
%       only equality constraints, no inequality constraints. Specifically,
%       this algorithm solves the optimization problem: 
%        minimize_x x'*B*x+2*a'*x
%        subject to A*x=b
%       If B is positive semidefinite, then if a solution exists, an
%       infinite number actually exist. This function returns one solution
%       and a matrix that can be used to find the other solutions.
%
%INPUTS: B An nXn real positive definite or positive semidefinite,
%          symmetric matrix.
%        a An nX1 real vector.
%        A A numConstXn real matrix. This can be an omitted or empty matrix
%          can be passed if there are no constraints. numConst<=n. There
%          shouldn't be any redundant constraints.
%        b A numConstX1 real vector.
%
%OUTPUTS: x The nX1 solution.
%         M An nXn matrix. All other solutions are of the form
%           xOther=x+M*y, where y is a an nX1 real vector.
%
%This function implements the equations of Chapter IV of [1]. Note that
%based on the constraints, it is possible for no solution to exist. This
%function is run assuming that a solution exists; it does not check for
%soluion non-existence.
%
%EXAMPLE:
%This is example 4.1 in [1]. x will be [9;18;35]/20.
% B=[1, 2,-1;
%    2, 4,-2;
%   -1,-2,1];
% a=[2;4;1];
% A=[1,2,1];
% b=4;
% [x,M]=semidefQuadProgEqConst(B,a,A,b)
%One sees that the third row and column of M is all zeros. Also,
%M(:,2)=-M(:,1)/2, so the matrix is singular. Thus, the second column of M
%can be used as a basis. So, to get any new  point, we can use x+M(:,2)*v
%(for a scalar v) or we could use x+5*M(:,2)*v, which is the same as the
%solution given in [1].
%
%REFERENCES:
%[1] D. L. Nelson, "Quadratic programming techniques using matrix
%    pseudoinverses," Ph.D. dissertation, Texas Tech University, Lubbock,
%    TX, Dec. 1969.
%
%October 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(A))
    %If there are no constraints, call the unconstrained algorithm.
    [x,M]=semidefQuadProgUnconst(B,a);
    return;
end

n=length(a);
I=eye(n,n);
pinvA=pinv(A);
diffA=(I-pinvA*A);
C=diffA*B*diffA;
pinvC=pinv(C);

%A combination of Equation 4.3 and 4.5.
x=pinvA*b+diffA*pinv(C)*diffA*(a-B*pinvA*b);
M=diffA*(I-pinvC*C);

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
