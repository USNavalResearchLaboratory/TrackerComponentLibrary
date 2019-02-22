function D=duplicationMatrix(n)
%%DUPLICATIONMATRIX Return a duplication matrix of a particular
%                   dimensionality. A duplication matrix D for the nXn
%                   symmetric matrix A is such that D*vech(A)=A(:). That
%                   is, the duplication matrix takes a vector of the
%                   diagonal and lower-triangular elements of A and returns
%                   A(:), essentially filling in the upper triangular part.
%
%INPUTS: n The number of dimensions of the square, symmetric matrix that is
%          to be duplicated.
%
%OUTPUTS: D An (n^2)X(n*(n+1)/2) duplication matrix. THis just consists of
%           zeros and 1s.
%
%Basic properties of the duplication matrix are given in [1].
%
%REFERENCES:
%[1] J. R. Magnus and H. Neudecker, "The elimination matrix: Some lemmas
%    and applications," SIAM Journal on Algebraic Discrete Methods, vol. 1,
%    no. 4, pp. 422-449, 1980.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

D=zeros(n^2,n*(n+1)/2);
idx=n+1;
row=0;

for i=1:n
    idx=idx-1;
    n2=n;
    
    col1=n-idx+1;
    if(idx~=n)
        for k=(idx+1):n
            row=row+1;
            D(row,col1)=1;
            n2=n2-1;
            col1=col1+n2;
        end
    end
    
    for j=1:idx
       row=row+1;
       
       col2=col1+j-1;
       D(row,col2)=1;
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
