function Rd=diagCovMatCoveringWorstDiags(rDiag)
%%DIAGCOVMATCOVERINGWORSTDIAGS Consider the ellipsoid defined by
% (z-zHat)'*inv(R)*(z-zHat)<=gamma
% However, suppose that only the rDiag=diag(R) values are known. We wish to
% create a diagonal matrix Rd such that (z-zHat)'*inv(Rd)*(z-zHat)<=gamma
% for all points z in the first ellipse. Since the cross terms are unknown,
% one must inflate the rDiag to guarantee that the origonal ellipsoid is
% covered regardless of the cross terms.
%
%INPUTS: rDiag A numDimx1 or 1XnumDim vector of the diagonal elements of
%              the matrix.
%
%OUTPUTS: Rd A numDimXnumDim diagonal matrix that has its diagonal scaled
%            to guarantee that it covers an ellipse with a matrix R having
%            diagonals rDiag, regardless of the cross terms.
%
%One can go through the steps that went into the derivation in
%diagCovMatCoveringMat. For the normalized matrix with 1's on the diagonal,
%the question is, what could the maximum eigenvalue be? The answer is
%numDim. This is the case where the matrix is singular, so the trace of the
%matrix equals the maximum eigenvalue. Thus, as diagCovMatCoveringMat
%explains how to inflate a covariance matrix, here, the answer turns out to
%just multiply diag(rDiag) by the number of dimensions.
%
%EXAMPLE:
%This plots many ellipses in 2D with a given set of diagonal elements and
%varying correlations. Then, we plot the ellipse for the original diagonal
%matrix, which is too small and finally, we plot the matrix from the output
%fo this function, which covers the worst-case correlation scenario.
% figure(1)
% clf
% hold on
% gammaVal=3;
% d=[2;1.1];%The diagonal elements.
% numPts=50;
% rho=linspace(-0.999,0.999,numPts);
% for k=1:numPts
%     diagVal=sqrt(d(1)*d(2))*rho(k);
%     R=[d(1),diagVal;
%        diagVal, d(2)];
%     drawEllipse([0;0],inv(R),gammaVal)
% end
% RDiag=diag(d);
% RFit=diagCovMatCoveringWorstDiags(d);
% drawEllipse([0;0],inv(RDiag),gammaVal,'linewidth',2)
% drawEllipse([0;0],inv(RFit),gammaVal,'linewidth',2)
% axis([-4,4,-4,4])
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=length(rDiag);
Rd=numDim*diag(rDiag);

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
