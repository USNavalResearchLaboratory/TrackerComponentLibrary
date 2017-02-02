function C=KronColumnwise(A,B)
%%KRONCOLUMNWISE Evaluate the Khatri-Rao product of two matrices taking the
%        columns as the partitions of the two matrices. Given a qXu matrix
%        A and a tXu matrix B, this product is the concatenation of
%        kron(A(:,i),B(:,i)) for i=1 to u. Thus, the resulting matrix has
%        dimensions .
%
%INPUTS: A A qXu matrix.
%        B A tXu matrix.
%
%OUTPUTS: C The (t*q)Xu Khatri-Rao product of A and B takign he partitions
%           column-wise.
%
%This product arises in a number of identities involving matrix calculus,
%such as in [1].
%
%REFERENCES:
%[1] J. W. Brewer, "Kronecker products and matrix calculus in systems
%    theory," IEEE Transactions on Circuits and Systems, vol. CAS-25, no.
%    9, pp. 772-781, Sep. 1978.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

q=size(A,1);
u=size(A,2);
t=size(B,1);

C=zeros(q*t,u);

for k=1:u
    C(:,k)=kron(A(:,k),B(:,k));
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
