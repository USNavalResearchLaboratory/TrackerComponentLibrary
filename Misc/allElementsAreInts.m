function boolVal=allElementsAreInts(x)
%%ALLELEMENTSAREINTS Given a vector or a matrix, this function returns true
%            if all of the elements in the matrix are exactly integers
%            (including complex integers. This differs from the isinteger
%            function that is built into Matlab as this actually checks the
%            values of the data, whereas isinteger only tells whether the
%            data type is an integer. THis function can handle floating
%            point numbers.
%
%INPUTS: x A scalar or matrix/hypermatrix where one wants to determine
%          whether all elements are integers.
%
%OUTPUTS: boolVal This is true if all elements of x are integer values
%               (regardless of data type) and this is false otherwise.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

boolVal=all(x(:)==fix(x(:)));

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
