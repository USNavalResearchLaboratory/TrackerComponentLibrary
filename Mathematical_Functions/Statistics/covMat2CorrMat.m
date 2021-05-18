function [C,sigma]=covMat2CorrMat(X)
%%COVMAT2CORRMAT Given a covariance matrix, compute the Pearson correlation
%                matrix. The essentially finds a matrix C with ones on the
%                main diagonal such that X=diag(sigma)*C*diag(sigma), where
%                sigma is the 
%
%INPUTS: X An nXn covariance matrix. This can be real and symmetric or
%          complex and Hermitian.
%
%OUTPUTS: C The nXn Pearson correlation matrix associated with X. The main
%           diagonal is all ones. 
%
%EXAMPLE:
%Here, we create a correlation matrix from a covariance matrix and then
%transform the correlation matrix back into a covariance matrix.
% X=[ 55,  -7, -12, -17;
%     -7,  45,  17,  22;
%    -12,  17,  35,  27;
%    -17,  22,  27,  25];
% [C,sigma]=covMat2CorrMat(X)
% XBack=diag(sigma)*C*diag(sigma)
%One will see that XBack is (within finite precision limitations) the same
%as X.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Get the vector of standard deviations.
sigma=sqrt(diag(X));

D=diag(1./sigma);
C=D*X*D;
%Force symmetry
C=(C+C')/2;
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
