function [QFactForm,R,Q]=HouseholderQR(A)
%%HOUSEHOLDERQR Perform a QR decomposition on A using a Householder
%           transformation. That is, find Q and R such that A=Q*R, where Q
%           is an orthogonal matrix and R is an upper-triangular matrix.
%           Unlike Matlab's standard QR function, this returns Q, though
%           the standard Q can also be obtained.
%
%INPUTS: A The mXn matrix whose QR decomposition is desired. A can be real
%          or complex.
%
%OUTPUTS: QFactForm The Q matrix of the decomposition in factorial form.
%           This is an mXn matrix. To evaluate Q*C, where C is a matrix,
%           one can use QRFactFormMult(QFactForm,C). Performing
%           multiplication this way rather than using Q (which is derived
%           from the matrix) is more efficient.
%         R The R matrix in the decomposition. If m>=n, then R is nXn. If
%           m<n, then R is mXn. For m>=n, the R matrix is essentially the
%           same as the "compact" QR decomposition of Matlab's QR function.
%           To make it non-compact, just appending m-n rows of zeros to the
%           bottom.
%         Q If requested, the actual Q matrix from the QR decomposition can
%           be obtained.
%
%Algorithm 5.2.1 of Chapter 5.2.2 of [1] is used. The function
%HouseholderVec also supports the complex Householder transformation,
%making this function valid for complex matrices A. To recover the actual Q
%matrix from the Q matrix in factorial form, Equation 5.1.5 of Chapter
%5.1.6 is used. The algorithms have been slightly modified so that they can
%be used with matrices where m<n (more columns than rows).
%
%If m<=n, then from the decomposition, one finds that A=Q*R and
%A=QRFactFormMult(QFactForm,R)
%If m>n, one still finds that A=Q*R, but this time,
%A=QRFactFormMult(QFactForm,[R;zeros(m-n,n)])
%as R is provided in a "compact" format.
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%November 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

m=size(A,1);
n=size(A,2);

for j=1:min(n,m-1)
    [v,beta]=HouseholderVec(A(j:m,j));
    A(j:m,j:n)=A(j:m,j:n)-(beta*v)*(v'*A(j:m,j:n));
    A((j+1):m,j)=v(2:(m-j+1));
end

%This comparison  and change of shape is necessary to handle matrices with
%m<n compared to the text.
if(m>=n)
    R=triu(A(1:n,1:n));
else
    R=triu(A(1:m,1:n));
end
QFactForm = tril(A,-1);

%Equation 5.1.5. If the non-factorial form of Q is desired. This has been
%modified from that of the paper to deal with the m<n case, in which case Q
%is not square.
if(nargout>2)
    if(m>=n)
        Q=eye(m,n);
    else
        Q=eye(m,m);
    end

    for j=min(n,m-1):-1:1
        v(j:m)=[1;QFactForm((j+1):m,j)];
        beta=2/(1+QFactForm((j+1):m,j)'*QFactForm((j+1):m,j));
        Q(j:m,j:end)=Q(j:m,j:end)-(beta*v(j:m))*(v(j:m)'*Q(j:m,j:end));
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
