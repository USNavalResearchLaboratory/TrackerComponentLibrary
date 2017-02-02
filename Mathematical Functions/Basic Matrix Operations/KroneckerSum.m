function kSum=KroneckerSum(A,B)
%KRONECKERSUM Find the Kronecker sum of an nXn matrix A and an mXm matrix
%             B. This is defined as kron(A,eye(m))+kron(eye(n),B)
%
%INPUTS: A An nXn matrix.
%        B An mXm matrix.
%
%OUTPUTS: kSum The (n*m)X(n*m) Kronecker sum of A and B.
%
%The Kronecker sum arises in numerous matrix identities. For example, it is
%used  in [1] for expressions for eigenvectors. 
%
%REFERENCES:
%[1] J. W. Brewer, "Kronecker products and matrix calculus in systems
%    theory," IEEE Transactions on Circuits and Systems, vol. CAS-25, no.
%    9, pp. 772-781, Sep. 1978.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

n=size(A,1);
m=size(B,1);

kSum=kron(A,eye(m))+kron(eye(n),B);

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
