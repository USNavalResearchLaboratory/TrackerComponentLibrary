function pinvX=pseudoInverse(X,algorithm)
%%PSEUDOINVERSE Compute the matrix pseudoinverse using a singular value
%               decomposition. This is essentially equivalent to the pinv
%               function in Matlab, except one can choose the algorithm to
%               use for determining the rank of the matrix.
%
%INPUTS:  X An NXM matrix whose pseudoinverse is desired.
% algorithm This is an optional input that specified the algorithm used in
%           matrixRank for determining the rank of the matrix. This affects
%           a threshold for declaring small singular values to be zero. The
%           default if omitted or an empty matrix is passed is 0.
%
%OUTPUTS: pinvX The MXN pseudoinverse of X.
%
%The definition of the pseudoinverse in terms of a singular value
%decomposition is given in Chapter 5.5.2 of [1].
%
%REFERENCES:
%[1] G. H. Golub and C. F. Van Loan, Matrix Computations, 4th ed.
%    Baltimore: Johns Hopkins University Press, 2013.
%
%June 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(algorithm))
    algorithm=0;
end

[rankVal,V,U,S]=matrixRank(X,algorithm);
%The singular values in S are in descending order.
s=diag(S);
s(1:rankVal)=1./s(1:rankVal);
s((rankVal+1):end)=0;

pinvX=V*diag(s)*U';

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
