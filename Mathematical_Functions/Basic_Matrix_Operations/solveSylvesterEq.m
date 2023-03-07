function X=solveSylvesterEq(A,B,C)
%%SOLVESYLVESTEREQ Solve the Sylester equation A*X+X*B=C for X, where A, B,
%           and C can be real or complex. This implementation only works
%           when a unique solution for X exists, which is the case when
%           none of the eigenvalues of A and B when summed together equals
%           0.
%
%INPUTS: A A pXp matrix.
%        B An rXr matrix.
%        C A pXr matrix.
%
%OUTPUTS: X The pXr matrix solving A*X+X*B=C.
%
%This function implements the Bartels-Stewart algorithm of [1]. The back
%substitution algorithm of Algorithm 7.6.2 in [2] is used for the second
%half of the algorithm, which also means that upper-triangular Schur
%decompositions are used for the first half of the algorithm.
%
%EXAMPLE:
%In this example, we solve a Sylvester equation and show that the result
%Tolerance is about 7e-14, which is good agreement.
% A=[3,  0,   0, 23;
%   -2,  5, -14, -6;
%   11, 11,  -7,  7;
%  -10, 15, -10, -1];
% B=[8,  4,   2,  6;
%   -7, -1,   1,  8;
%  -14, -1,  15, -2;
%  -14, 14,  -8,  2];
% C=[16,  3,  6,  8;
%     3, 11,  8, 11;
%     6,  8,  6, 13;
%     8, 11, 13,  1];
% X=solveSylvesterEq(A,B,C);
% AbsTol=max(max(abs(A*X+X*B-C)))
%
%REFERENCES
%[1] R. H. Bartels and G. W. Stewart, "Algorithm 432: Solution of the
%    matrix equation AX+XB=C [F4]," Communications of the ACM, pp. 820-826,
%    Sep. 1972.
%[2] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: The Johns Hopkins Press, 2013.
%
%December 2022 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%We use the complex Schur decomposition. If we used the real one, then
%Ap and Bp are not guaranteed to be upper triangular if the matrices have
%complex eigenvalues.
[U,Ap]=schur(A,'complex');
[V,Bp]=schur(B,'complex'); 
Cp=U'*C*V;

%Solve the transformed system using Algorithm 6.7.2 of [1].
Bp=-Bp;

r=size(Ap,1);
p=size(Bp,1);

I=eye(r,r);

for k=1:r
    Cp(1:p,k)=Cp(1:p,k)+Cp(1:p,1:(k-1))*Bp(1:(k-1),k);
    z=(Ap-Bp(k,k)*I)\Cp(1:p,k);
    Cp(1:p,k)=z;
end

%Cp hold the transformed result. Now undo the transformation.
X=U*Cp*V';

if(isreal(A)&&isreal(B)&&isreal(C))
    %If all of the inputs are real, then Cp should be real, not accounting
    %for finite precision limits.
    X=real(X);
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
