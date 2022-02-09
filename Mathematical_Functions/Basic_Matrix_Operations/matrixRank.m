function [rankVal,V,U,S]=matrixRank(A,algorithm)
%%MATRIXRANK Determine the rank of a matrix using a chosen criterion.
%
%INPUTS: A The MXN matrix whose rank is desired.
% algorithm An optional parameter specifying the criterion to use. All
%          methods compare the singular values to a bound. Possible values
%          are:
%          0 (The default if omitted or an empty matrix is passed) Use
%            eps(norm(A,1)) as the bound.
%          1 Use max(size(A))*eps(max(s)), where s is the vector of
%            singular values.
%          2 Use eps()*norm(A,1).
%          3 Use max(size(A))*eps()*max(s).
%
%OUTPUTS: rankVal The rank of matrix A.
%               V The right singular vectors from the SVD used here (The
%                 full SVD is computed, not the 'econ' or 0 options). Given
%                 the rank, one can extract a basis or a nullspace.
%               U The left singular vectors from the SVD.
%               S A diagonal matrix of singular values from the SVD.
%
%The notion of the matrix rank relating to a thresholding of the singular
%values is discussed in Chapter 5.4.1 of [1].
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%September 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

if(nargout>1)
    [U,S,V]=svd(A);
else
    [~,S,~]=svd(A);
end
s=diag(S);

switch(algorithm)
    case 0
        rankVal=sum(s>eps(norm(A,1)));
    case 1
        rankVal=sum(s>max(size(A))*eps(max(s)));
    case 2
        rankVal=sum(s>eps()*norm(A,1));
    case 3
        rankVal=sum(s>max(size(A))*eps()*max(s));
    otherwise
        error('Unknown rank algorithm specified.')
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
