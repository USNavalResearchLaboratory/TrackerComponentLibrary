function RTilde=diagCovMatCoveringMat(R)
%%DIAGCOVMATCOVERINGMAT Consider a ellipsoidal region specified by
%           (z-zHat)'*inv(R)*(z-zHat)<=gamma, where gamma is a positive
%           threshold and R is not a diagonal matrix. Determine a diagonal
%           matrix Rd such that the ellipsoid
%           (z-zHat)'*inv(Rd)*(z-zHat)<=gamma completely contains the
%           original ellipsoid and also touches the original ellipsoid (but
%           does not intersect it).
%
%INPUTS: R An nXn real, positive (semi)definite symmetric matrix.
%
%OUTPUTS: RTilde A diagonal matrix specifying a probability region that
%            encompasses that specified by R.
%
%This finds a diagonal matrix to transform R to have 1's on the diagonal.
%Then, the maximum eigenvalue is found, which is related to the maximum
%distance of a point on the ellipse from the origin. A circle (diagonal
%matrix with all equal values) with the radius proportional to the maximum
%eigenvalue will cover the ellipse implied by R. Thus, the original scaling
%is undone so that it still works in the original coordinate system.
% 
%EXAMPLE:
%This is an example in 2D. The original ellipse is shown along with what
%one would get by just dropping the cross terms as well as the output of
%this function. one can see that just dropping the cross terms leads to a
%matrix that is too small. However, the output of this function, completely
%covers the original ellipse. 
% R=[2,-1.3;
% -1.3, 1.1];
% RDiag=diag(diag(R));
% Rd=diagCovMatCoveringMat(R);
% gamma=16;
% figure(1)
% clf
% hold on
% drawEllipse([0;0],inv(R),gamma,'linewidth',2)
% drawEllipse([0;0],inv(RDiag),gamma,'linewidth',2)
% drawEllipse([0;0],inv(Rd),gamma,'linewidth',2)
% legend('With Cross Terms','Diagonal Only','Diagonal, Covering')
% axis([-10, 10, -10, 10])
%
%February 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

rDiag=diag(R);
DInv=diag(1./sqrt(rDiag));
lambdaMax=max(eig(DInv*R*DInv));
RTilde=lambdaMax*diag(rDiag);

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
