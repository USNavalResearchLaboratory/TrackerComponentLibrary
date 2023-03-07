function [d,F,G]=principalAngVec(A,B)
%%PRINCANGPRINCVEC Given an mXp matrix and an mXq matrix with p>=q, compute
%           the principal angles and vectors between the matrices. This can
%           be used in canonical correlation analysis to find combinations
%           of variables with maximum correlation between each other.
%           Specifically, the principle angles are defined recursively by
%           d(k)=cos(theta_k)=F(:,k)'*G(:,k)=
%               max_{F(:,k)} max_{G(:,k)} F(:,k)'*G(:,k)
%           where in the optimization F(:,k)  is a unit vector in the space
%           spanned by A that is orthogonal to all F(:,i) with i<k, and
%           G(:,k) is similarly a unit vector in the space spanned by B
%           that is orthogonal to all G(,i) with i<k.
%
%INPUTS: A An mXp matrix with independent columns.
%        B An mXq matrix with independent columns q<=p.
%
%OUTPUTS: d A qX1 vector of the cosine values of the principle angles (Use
%           acos to get the angles). The angles all satisfy 0<=theta<=pi/2.
%         F mXq The first set of principal vectors.
%         G mXq The second set of principal vectors.
%
%This function implements Algorithm 6.4.3 in Chapter 6.4 of [1].
%
%EXAMPLE:
%This function can be used to quickly tell if two matrices are identicle
%except for a reordering and scaling of the columns. If they are, then the
%principal angles will all be zero (the cosines will all be 1). We
%demonstrate this here: 
% A=[8,0,1;
%    9,2,9;
%    1,5,9;
%    9,9,4;
%    6,9,8];
% B=-4*A(:,[3,1,2]);%Scale and permute the columns.
% [d,F,G]=principalAngVec(A,B)
%One will see that the three d values are all 1 within finite precision
%bounds.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

q=size(B,2);

if(q>size(A,2))
    error('The number of columns of B must be <= the number of columns of A.') 
end

[QA,~]=qr(A,0);
[QB,~]=qr(B,0);
[Y,D,Z]=svd(QA'*QB);
F=QA*Y(:,1:q);
G=QB*Z(:,1:q);
d=diag(D);

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
