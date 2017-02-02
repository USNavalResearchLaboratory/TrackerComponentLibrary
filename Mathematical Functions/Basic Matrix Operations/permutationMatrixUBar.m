function UBar=permutationMatrixUBar(p,q)
%%PERMUTATIONMATRIXUBAR Obtain the permutation matrix UBar. This is a
%        specific type of permutation matrix that plays a role in matrix
%        calculus. The partial derivative of a pXq matrix A with respect to
%        itself equals permutationMatrixUBar(p,q) as in [1]. This matrix is
%        defined as 
%        UBar=\sum_{i=1}^p\sum_{k=1}^q kron(Epq(i,k),Epq(i,k))
%        where Epq(i,k) returns a pXq matrix of zeros with a 1 in row i and
%        column j.
%
%INPUTS: p,q The dimensions defining the ultimate size of the resulting
%            permutation UBar matrix.
%
%OUTPUTS: UBar A p*q permutation matrix having the above mentioned
%              properties.
%
%This function just implements Equation 5 in [1].
%
%This type of a permutation matrix plays a notable role in matrix calculus
%used in control systems [1].
%
%REFERENCES:
%[1] J. W. Brewer, "Kronecker products and matrix calculus in systems
%    theory," IEEE Transactions on Circuits and Systems, vol. CAS-25, no.
%    9, pp. 772-781, Sep. 1978.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

UBar=zeros(p^2,q^2);

Epq=zeros(p,q);
for i=1:p
    for k=1:q
        Epq(i,k)=1;
        
        UBar=UBar+kron(Epq,Epq);
        
        Epq(i,k)=0;
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
