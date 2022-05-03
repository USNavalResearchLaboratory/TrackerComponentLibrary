function Q=factFormMat2Explicit(A)
%%FACTFORMMAT2EXPLICIT Sequences of Householder transformations are often
%           stored in a factored form rather than an explicit form. This
%           function takes the factored form and converts it to the
%           explicit form. As described in Chapter 5.1.6 of [1], the
%           explicit form of the transformation is given by the product of
%           matrices: Q=A1*A2*A3...*An where each Ai matrix is of the form
%           Ai=eye(m)-betai*(vi*vi') and the matrix A is of size mXn. The
%           scalar beta value can be recovered from the v vector (see
%           below), so the vector vi encodes all of the information in the
%           matrix product. The vector vi is of the form
%           [zeros(i-1);1;vFull(m-i)]. The ith column of the matrix A holds
%           the values of vFull(m-i), which is a length m-i vector. The
%           value are put in the lower triangular portion of A, not
%           including the diagonal. such vectors are often obtained using
%           the HouseholderVec function.
%
%INPUTS: A A mXn matrix formatted as described above.
%
%OUTPUTS: Q The explicit form of the product of subsequent Householder
%           matrices A1*A2*...
%
%Chapter 5.1.6 describes the factored form representation. It often arises
%in QR and bidiagonalization algorithms. Equation 5.1.5 gives the
%conversion. However, it had to be slightly modified to deal with a the
%vase when a particular column of A is all zeros.
%
%For the instance when vFull equals zero, that is when the vi vector above
%is of the form [zeros(i-1);1;zeros(m-i)], one can look to Algorithm 5.1.1
%in Chapter 5.1.3 of [1] to see that occurs either when beta=0 or beta=2
%(correcting using the errata for the book, which listed -2). However, one
%doesn't know which branch to take. Here, we assume that the factored form
%was computed using something like HouseholderVec with forceSign=false.
%This is because Algorithm 5.1.3 does NOT account for the possible beta=0
%branch in Algorithm 5.1.1 and just from the factored form, there is no way
%to know which branch was taken.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%January 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(A,1);
n=size(A,2);

%This if-statement deals with the case where m<n, which is not explicitly
%given in the book.
if(m>=n)
    Q=eye(m,n);
else
    Q=eye(m,m);
end
    
v=zeros(m,1);
%The min operation is not explicitly given in the book. the m-1 deals with
%j+1 being past the end of the rows of the matrix and the n as a minimum
%deals with j not being past the end of the columns of the matrix.
for j=min(n,m-1):-1:1
    AVec=A((j+1):m,j);

    v(j:m)=[1;AVec];
    normAVec2=AVec'*AVec;
    beta=2/(1+normAVec2);
    Q(j:m,j:end)=Q(j:m,j:end)-(beta*v(j:m))*(v(j:m)'*Q(j:m,j:end));
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
