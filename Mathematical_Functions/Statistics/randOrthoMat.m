function randOrthoMat=randOrthoMat(numDim)
%RANDORTHOMAT Generate a random orthogonal matrix of the desired
%             dimensionality (with respect to the Haar distribution over
%             the orthogonal group --kind of a "uniform" distribution).
%             This means that the dot product of any two columns is zero.
%             The orthogonal matrices generated have determinants of +1 or
%             -1.
%
%INPUTS: numDim The scalar number of dimensions that the random orthogonal
%               matrix should have.
%
%OUTPUTS: randOrthoMat A numDimXnumDim random orthogonal matrix with a
%                     determinant of either +1 or -1.
%
%As noted in [1], random orthogonal matrices can be obtained by first
%gennerating a matrix of independent standard normal random variables, then
%takign the QR decomposition of the matrix such that R has all positive
%elements along the diagonal. However, since Matlab's QR command can
%produce R matrices with negative diagonal elements, we have to effectively
%transfer the sign of R onto Q to make the elements of the diagonal of R
%all positive. This means multiplying Q by diag(sign(diag(R))).
%
%EXAMPLE:
%We will plot the random directions in 2D on the unit circle to show that
%this produces a relatively uniform set of rotation matrices.
% numPoints=500;
% figure(1)
% clf
% hold on
% for curPoint=1:numPoints
%     thePoint=randOrthoMat(2)*[1;0];
%     scatter(thePoint(1),thePoint(2),'.k','linewidth',2)
% end
%
%REFERENCE:
%[1] G. W. Stewart, "The efficient generation of random orthogonal matrices
%    with an application to condition estimators," SIAM Journal on
%    Numerical Analysis, vol. 17, no. 3, pp. 403-409, Jun. 1980.
%
%October 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

X=randn(numDim,numDim);
[Q,R]=qr(X);

randOrthoMat=Q*diag(sign(diag(R)));

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
