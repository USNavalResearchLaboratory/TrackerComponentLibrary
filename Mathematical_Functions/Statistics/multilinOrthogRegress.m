function [n,a]=multilinOrthogRegress(zSamp)
%%MULTILINORTHOGREGRESS Perform multilinear regression using orthogonal
%            distances as the penalties. This is also known as Deming
%            regression. It is a type of total least squares estimation.
%            Here, the equation for a line/plane/hyperplane/etc. is taken
%            to be dot(n,x-a)=0, where n is a unit vector and x is the free
%            parameter. This function finds n and a to minimize
%            sum_{i=1}^numPoints dot(n,zSamp(:,i)-a).
%
%INPUTS: zSamp A numDimsXnumPoints set of numPoints data samples.
%              numDims>=2.
%
%OUTPUTS: n The fitted numDimsX1 unit vector in the linear equation above.
%         a The fitted constant in the equation above.
%
%Formulae for this type of regression are given in Appendix A.7.4 of [1]
%and are implemented here.
%
%EXAMPLE:
% n=[1;2;3]/norm([1;2;3]);
% a=[-3;6;8];
% %The equation is dot(n,(z-a))=0. We want some z values. We will generate
% %some values in the first two dimensions on a grid and then get the third
% %dimension by solving the linear equation to get
% % z=(-n(1).*x-n(2).*y+dot(n,a))/n(3)
% numPts=10;
% xPts=linspace(-4,4,numPts);
% yPts=linspace(-4,4,numPts);
% [X,Y]=ndgrid(xPts,yPts);
% Z=(-n(1).*X(:)-n(2)*Y(:)+dot(n,a))/n(3);
% zTrue=[X(:).';Y(:).';Z(:).'];
% [nEst,aEst]=multilinOrthogRegress(zTrue)
% %One will see that nEst equals n within finite precision bounds. However,
% %aEst will not equal a. This is because a is not unique. However, one can
% %verify that dot(n,zTrue(:,i)-a) is zero, within finite precision bounds
% %for all i.
% zSamp=zTrue+0.01*randn(size(zTrue));
% [nEst1,aEst1]=multilinOrthogRegress(zSamp)
% %With noise added to both components, one gets a solution that is close
% %to the correct solution.
%
%REFERENCES:
%[1] P. J. Schneider and D. H. Eberly, Geometric Tools for Computer
%    Graphics. Amsterdam: Morgan Kaufmann Publishers, 2003.
%
%April 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    a=mean(zSamp,2);

    y=bsxfun(@minus,zSamp,a);
    M=y*y';
    [V,d]=eig(M,'vector');
    [~,minIdx]=min(d);
    n=V(:,minIdx);
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
