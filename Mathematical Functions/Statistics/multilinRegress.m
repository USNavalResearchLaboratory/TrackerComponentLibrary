function coeffs=multilinRegress(zSamp)
%%MULTILINREGRESS Perform multilinear regression. This means, given a set
%          of multidimensional points of the form [x1;x2;x3;...;xn],
%          determine the line/plane/hyperplane best fitting the points
%          such that
%          xn=coeffs(1)*x2+coeffs(2)*x2+....coeffs(n-1)*x(n-1)+coeffs(n).
%          The "best" fitting line/plane/hyperplane is the one that
%          minimizes
%          sum_{i=1}^n (dot(coeffs(1:(n-1)),zSamp(1:(n-1),i))+coeffs(n)-zSamp(n))^2
%
%INPUTS: zSamp A numDimsXnumPoints set of numPoints data samples.
%              numDims>=2.
%
%OUTPUTS: coeffs The numDimsX1 set of weights for the linear equation
%                described above.
%
%Formulae for this type of regression are given in Appendix A.7.1 and A.7.3
%of [1] for 2D and 3D measurements. This function just continues the
%pattern implementing the results for an arbitrary number of dimensions.
%Note that this function will fail if the plane being estimated does not
%depend on the final parameter. In 2D, this would correspond to a vertical
%line (e.g. x=4).
%
%EXAMPLE:
% coeffs=[4;6;2];
% numPts=20;
% xVals=linspace(-4,4,numPts);
% yVals=linspace(-4,4,numPts);
% [x,y]=ndgrid(xVals,yVals);
% xyPts=[x(:).';y(:).'];
% zPts=sum(bsxfun(@times,coeffs(1:2),xyPts),1)+coeffs(3);
% zPts=zPts+0.1*randn(1,numPts^2); %Add noise.
% zSamp=[xyPts;zPts];
% coeffsEst=multilinRegress(zSamp)
%One finds that coeffsEst is close to coeffs. A perfect fit is achieved if
%there is zero noise.
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDims=size(zSamp,1);
numPoints=size(zSamp,2);

temp=[zSamp(1:(numDims-1),:);ones(1,numPoints)];
A=temp*temp';
b=sum([bsxfun(@times,zSamp(1:(numDims-1),:),zSamp(numDims,:));zSamp(numDims,:)],2);

coeffs=A\b;
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
