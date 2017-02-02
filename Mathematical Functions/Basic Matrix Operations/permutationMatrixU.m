function U=permutationMatrixU(p,q)
%%PERMUTATIONMATRIXU Obtain the permutation matrix U. This is a specific
%        type of permutation matrix that plays a role in matrix calculus
%        and in equations involving Kronecker products. Specifically,
%        U=\sum_{i=1}^p\sum_{k=1}^q kron(Epq(i,k),Eqp(k,i))
%        where Epq(i,k) returns a pXq matrix of zeros with a 1 in row i and
%        column j and Eqp(k,i) similarly returns a qXp matrix with a 1 in
%        row k and column i.This permutation matrix also arises when taking
%        matrix derivatives. Table V of [1] shows that the partial matrix
%        derivative of the transpose of pXq matrix A with respect to itself
%        is permutationMatrixU(p,q), assuming all elements of q are
%        independent.
%
%INPUTS: p,q The dimensions defining the ultimate size of the resulting
%            permutation U matrix.
%
%OUTPUTS: U A (p*q)X(p*q) permutation matrix having the above mentioned
%           properties.
%
%This function just implements Equation 4 in [1]. A major application of
%this matrix is based on the identity that
%kron(B,A)=U1*kron(A,B)*U2 where U1=permutationMatrixU(s,p),
%U2=permutationMatrixU(q,t) and A and B are respectively pXq and sXt
%matrices.
%
%This type of a permutation matrix plays a notable role in matrix calculus
%used in control systems [1] and a matrix representation of general
%relativity theory [2]. Specifically, in [2], the matrix
%permutationMatrixU(4,4) is often used.
%
%REFERENCES:
%[1] J. W. Brewer, "Kronecker products and matrix calculus in systems
%    theory," IEEE Transactions on Circuits and Systems, vol. CAS-25, no.
%    9, pp. 772-781, Sep. 1978.
%[2] G. Ludyk, Einstein in Matrix Form: Exact Derivation of the Theory of
%    Special and General Relativity without Tensors. Heidelberg: Springer,
%    2013.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

U=zeros(p*q,p*q);

Epq=zeros(p,q);
Eqp=zeros(q,p);
for i=1:p
    for k=1:q
        Epq(i,k)=1;
        Eqp(k,i)=1;
        
        U=U+kron(Epq,Eqp);
        
        Epq(i,k)=0;
        Eqp(k,i)=0;
    end
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
