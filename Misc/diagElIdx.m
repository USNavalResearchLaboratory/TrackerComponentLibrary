function idx=diagElIdx(N,minorDiag)
%%DIAGELIDX For an NXN matrix, determine the linear indices of the elements
%           on the main or minor diagonal of the matrix.
%
%INPUTS:  N The size of the NXN matrix.
% minorDiag Optionally, this can specify whether the main or minor diagonal
%           of the matrix is desired. minorDiag=false (the default if
%           omitted or an empty matrix is passed) chooses the main
%           diagonal, and minorDiag=true chooses the minor diagonal of the
%           matrix.
%
%OUTPUTS: idx A 1XN vector of the indices on the main (or minor) diagonal
%             of the matrix. Minor diagonal elements go from the lower-left
%             corner of the matrix to the upper-right corner.
%
%The major diagonal is also known as the principal diagonal, the primary
%diagonal, the major diagonal and the leading diagonal. The minor diagonal
%of a matrix is also known as the antidiagonal, the counterdiagonal, the
%trailing diagonal, and the secondary diagonal.
%
%EXAMPLE:
% N=8;
% M=magic(N)
% mainDiagEls=M(diagElIdx(N))
% minorDiagEls=M(diagElIdx(N,true))
%One will see that mainDiagEls contains the main diagonal elements of M,
%and minorDiagEls contains the minor diagonal elements, which increase
%sequentialy.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(minorDiag))
    minorDiag=false;
end

if(minorDiag==false)
    idx=1:(N+1):(N^2);
else
    idx=N:(N-1):(N^2-N+1);
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
