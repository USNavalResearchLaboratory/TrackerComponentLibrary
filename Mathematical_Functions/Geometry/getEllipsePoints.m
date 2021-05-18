function zPoints=getEllipsePoints(z,A,gammaVal,numPoints,invertA)
%%GETELLIPSEPOINTS Generate numPoints around a 2D ellipse. These can be
%           used, for example, to plot the ellipse. A point zp on the
%           ellipse/ ellipsoid satisfies the equation
%           (zp-z)'*A*(zp-z)=gammaVal
%           or is invertA is true, the equation 
%           (zp-z)'*inv(A)*(zp-z)=gammaVal
%           gammaVal determines what amount of probability is in the
%           region. If omitted, gammaVal=16.2235, which corresponds
%           approximately to the 99.97% confidence region, is used (The
%           value comes from ChiSquareD.invCDF(0.9997,2)).
%
%INPUTS: z A 2XN vector corresponding to the centers of the N ellipses for
%          which points should be obtained.
%        A A 2X2XN set of N positive definite matrices that specify the
%          size and shape of the ellipses, where a point zp is on the ith
%          ellipse if
%          (zp-z(:,i))'*A(:,:,i)*(zp-z(:,i))=gammaVal (if invertA is true,
%          then replace A(:,:,i) with inv(A(:,:,i))).
% gammaVal An optional parameter specifying the size of the ellipse/
%          ellipsoid. If omitted or an empty matrix is passed, then
%          gammaVal=16.2235 is used. This is approximately the value for a
%          99.97% confidence region if A are inverse covariance matrices of
%          a Gaussian distribution. gammaVal must be positive. 
% numPoints An optional parameter specifying how many points should be
%          generated. The default if omitted or an empty matrix is passed
%          is 2000. If this number is not an even number, then the next
%          largest even number of points is generated.
%  invertA If this is true, then A is inverted before use. The default if
%          omitted or an empty matrix is passed is false.
%
%OUTPUTS: zPoints A 2XnumPointsXN set of the points for each of the
%                 ellipses.
%
%An eigendecomposition of a matrix A breaks it into parts V and D such that
%V*D*V'=A, where V is a rotation matrix and D is a diagonal matrix.
%Considering the problem here, this means that (zp-z)'*A*(zp-z)=
%(V'*zp-V'*z)'*D*(V'*zp-V'*z). If we substitute kp=V'*zp and k=V'*z, then
%the equation is (kp-k)'*D*(kp-k)=gammaVal, which is the equation for an
%ellipse where the axes are aligned with the coordinate axes. Thus, one can
%first find the points for an ellipse that is aligned with the coordinate
%axes (the kp points), and then rotate it back to the proper alignment (the
%zp points). However, it is simpler to first find an ellipse centered about
%zero. Thus, substitute l=kp-k and find the l points, then shift and rotate
%them to get the z-points.
%
%When considering an ellipse, suppose that l=[x;y] and D=diag(a,b). Then
%the ellipse equation is
%x^2*a+y^2*b=gammaVal
%Solving for y in terms of x, we have
%y=+/-sqrt((gammaVal-x^2*a)/b).
%The x values are limited so that the argument of the square root is
%positive. Thus the x values range from -sqrt(gammaVal/a) to
%+sqrt(gammaVal/a). Thus, a simple way to plot this centered, axis-aligned
%ellipse is to generate x-values in the valid range and then find the two
%sets of y values. This function does that and transforms the results to
%get the final result.
%
%EXAMPLE:
%Here, we generate points for two ellipses and then plot them.
% z=[[4;3],[3;4]];
% A=zeros(2,2,2);
% A(:,:,1)=[4,-2;
%           -2,5];
% A(:,:,2)=[5,-2;
%           -2,5];
% zPoints=getEllipsePoints(z,A,[],[],true);
% figure(1)
% clf
% hold on
% plot(zPoints(1,:,1),zPoints(2,:,1),'-r')
% plot(zPoints(1,:,2),zPoints(2,:,2),'-b')
%
%December 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(invertA))
    invertA=false;
end

if(nargin<4||isempty(numPoints))
    numPoints=2000; 
end

if(nargin<3||isempty(gammaVal))
    gammaVal=16.2235;
end 

if(invertA)
    A=applyFunToEachMatrix(@inv,A);
end

N=size(z,2);

numPoints=ceil(numPoints/2);

zPoints=zeros(2,2*numPoints,N);
for curEllips=1:N
    %Perform an eigenvalue decomposition of A.
    [V,D]=eig(A(:,:,curEllips));
    %We will put the smallest eigenvalue at the start, if it is
    %not already there.
    [~,idx]=sort(diag(D),'ascend');
    V=V(idx,idx);
    D=D(idx,idx);
    a=D(1,1);
    b=D(2,2);

    xBound=sqrt(gammaVal/a);
    x=linspace(-xBound,xBound,numPoints);
    %The real command deals with possible finite precision issues.
    y=real(sqrt((gammaVal-x.^2*a)/b));
    %The centered, axis-aligned ellipse points.
    l=[x,fliplr(x);
       y,-fliplr(y)];
    zp=bsxfun(@plus,V*l,z(idx,curEllips));
    zPoints(:,:,curEllips)=zp(idx,:);%Undo the reordering.
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
