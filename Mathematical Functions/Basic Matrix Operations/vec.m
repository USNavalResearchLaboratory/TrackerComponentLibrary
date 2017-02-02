function val=vec(A)
%%VEC Given A matrix, return the values of the matrix stacked columnwise.
%     This is just the same as A(:), but this function can be more
%     convenient to use when programming equations as one can
%     say things like c=b+vec(A*B) rather than temp=A*B; c=b+temp(:).
%
%INPUTS: A A matrix.
%
%OUTPUTS: val The values in the matrix stacked column-wise. If A is a
%             hypermatrix, then the values are stacked in order of the
%             indices of A (rows, columns, hypercolumns, etc.).
%
%February 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

val=A(:);
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
