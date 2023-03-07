function [P,r,c]=makeDoublyStochastic(A,tol,g,omitCheck)
%%DOUBLESTOCHBAL Transform a matrix into a doubly stochastic matrix. That
%                is, a matrix where all of the rows sum to one and all of
%                the columns sum to one. This function can be used to
%                battle precision errors that make doubly stochastic
%                matrices non-doubly stochastic after a number of
%                operations.
%
%INPUT:
% A: An NXN square matrix that has nonnegative entries and at least one
%    strictly positive diagonal as defined in [2].
% tol: An optional input that determines when the algorithm has converged.
%      This is a measure of the column sum errors. If tol is not provided
%      or is empty, the default value of eps*N*N is used.
% g: A small perturbation value that assures convergence of the algorithm
%    is the matrix A is of a type that is ill-suited to being balanced.
%    The default value if this parameter is omitted or empty is eps. If it
%    is known that the matrix A has the proper structure (as discussed in
%    the paper cited below), then eps can be set to zero.
% omitCheck: An optional boolean indicating whether to skip checking for
%            the positive diagonal condition. The default is false. Note
%            that if the check is skipped, the function may not reach a
%            stopping condition.
%
%OUTPUTS:
% P: A doubly stochastic matrix obtained by transforming A as
%    P=diag(r)*A*diag(c).
% r: One of the two vectors used to transform A into P.
% c: The second of the two vectors used to transform A into P.
%
%This function is an implementation of the algorithm given in [1], and the
%specific meaning of the values tol and g are discussed in the paper. The
%underlying Sinkhorn-Knopp algorithm was originally derived in [2].
%
%The implementation of the termination criterion of the paper has been
%slightly changed. The termination criterion given in the paper can not be
%used when entering the while-loop given in the paper to update c and r.
%For example, with tol=8.8818e-16 and A=[0.2,0.6;0.8,0.4], the termination
%criterion is fulfilled without entering the loop even though A is not
%doubly stochastic. Thus, the while-loop of the paper has been replaced
%with a do-while loop.

%A diagonal of matrix A as defined in [2] is a set containing the elements
%{A_(1,f(1)), A_(2,f(2)), ..., A_(1,f(m))} for a chosen permutation (shift)
%of indices {1,2,...,m}.
%
%EXAMPLE 1: Convert a magic square into a doubly-stochastic matrix. The
%           positive diagonal check can be omitted since all entries are
%           greater than 0. The row and column sums should be 0 to MATLAB's
%           default machine precision.
% Ain = magic(10);
% [Aout,d1,d2] = makeDoublyStochastic(Ain);
% assert(all(abs(1-sum(Aout,1))<tol))
% assert(all(abs(1-sum(Aout,2))<tol))
%
%EXAMPLE 2: Convert a random matrix into a doubly-stochastic matrix. The
%           seed given below will produce a matrix satisfying the
%           conditions of the Sinkhorn-Knopps algorithm. There will be a
%           noticeable deviation from 1 in the row and column sums, but
%           they will fall within the demanded tolerance.
% rng(42)
% Ain = rand(10);
% tol = 0.001;
% [Aout,d1,d2] = makeDoublyStochastic(Ain,tol);
% assert(all(abs(1-sum(Aout,1))<tol))
% assert(all(abs(1-sum(Aout,2))<tol))
%
%REFERENCES:
%[1] P. A. Knight, "The Sinkhorn-Knopp algorithm: Convergence and
%    applications," SIAM Journal on Matrix Analysis and Applications, vol.
%    30, no. 1, 2008.
%[2] R. Sinkhorn and P. Knopp, "Concerning nonnegative matrices and doubly
%    stochastic matrices," Pacific Journal of Mathematics, vol. 21, no. 2,
%    pp. 343-348, 1967.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%October 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(A,1);

if(nargin<2||isempty(tol))
    tol=eps*N*N;
end
if(nargin<3||isempty(g))
    g=eps;
end
if(nargin<4)
    omitCheck=false;
end

if N~=size(A,2)
    error('A must be a square matrix.')
end
if any(A<0)
    error('A must be a nonnegative matrix.')
end
if tol<0
    error('tol must be nonnegative.')
end

% Check for strictly positive diagonal to ensure solution exists.
hasPosDiag = false;
if nargin==3 && omitCheck
    hasPosDiag = true;
end
k = 0;
while ~hasPosDiag && k<N
    d = extendedDiag(A,k);
    if all(d>0)
        hasPosDiag = true;
    end
    k = k+1;
end
if hasPosDiag==false
    error('Ain must have at least one positive diagonal to ensure convergence.')
end

r=ones(N,1);
c=r;
d=A'*r+g*sum(r);
while(norm(c.*d-1,1)>tol)
   c=1./d;
   r=1./(A*c+g*sum(c));
   d=A'*r+g*sum(r);
end

P=diag(r)*A*diag(c);

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
