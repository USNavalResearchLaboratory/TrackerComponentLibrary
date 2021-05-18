function xyPoints=getEllipseHullPoints(z,A,gammaVal,numPoints,invertA)
%%GETELLIPSHULLPOINTS  Consider an ellipsoid in 3D (x,y,z) coordinates. For
%       every fixed z coordinate, if the x-y plane at the z coordinate
%       intersects the ellipsoid, it cuts an ellipse. For example, the
%       function projEllipse2ZPlane gets the parameters of the ellipse of
%       intersection. This function gets points marking the outer limits of
%       all of those ellipses in the x-y plane (the hull of the ellipses).
%
%INPUTS: z A 3XN vector corresponding to the centers of the N ellipses for
%          which points should be obtained.
%        A A 3X3XN set of N positive definite matrices that specify the
%          size and shape of the ellipsoids, where a point zp is on the ith
%          ellipsoid if
%          (zp-z(:,i))'*A(:,:,i)*(zp-z(:,i))=gammaVal (if invertA is true,
%          then replace A(:,:,i) with inv(A(:,:,i))).
% gammaVal An optional parameter specifying the size of the ellipse/
%          ellipsoid. If omitted or an empty matrix is passed, then
%          gammaVal=18.8049 is used. This is approximately the value for a
%          99.97% confidence region if A are inverse covariance matrices of
%          a Gaussian distribution. gammaVal must be positive gammaVal must
%          be positive.
% numPoints An optional parameter specifying how many points should be
%          generated. The default if omitted or an empty matrix is passed
%          is 500.
%  invertA If this is true, then A is inverted before use. The default if
%          omitted or an empty matrix is passed is false.
%
%OUTPUTS: xyPoints A 2XnumPointsXN set of the points for each of the
%                  ellipsoid projection hulls.
%
%The 3D ellipsoid (with the coordinate system shited to put the origin at
%the center) is defined as x'*R*x=gammaVal. We can rewrite this as 
%[r*u,z]*R*[r*u;z]=gammaVal
%where r is a positive scalar and u is a 2X1 unit vector. The matrix R can
%be broken up accordingly as
%R=[Rxy, rz;
%   rz', rzz];
%We want the maximum extent of the ellipsoid in the x-y plane for each
%direction (as specified by u). This will provide the hull of the ellipsoid
%points. Thus, we perform the optimization:
% maximize (over r,z) r^2
% such that r^2*u'*Rxy*u+2*r*z*rz'*u+z^2*rzz=gammaVal
%where the constraint is just from rewriting the expression for the
%definition of an ellipsoid. Writing the Lagrangian, (with lambda as the
%Lagrangian parameter) taking derivatives and setting them equal to 0, we
%get two additional equations. Solving those two equations with the
%original constraint, we get expressions for r, z, and lambda. Choosing the
%solution such that r is positive and nonzero, one gets the expressions
%implemented in this function.
%
%December 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(invertA))
    invertA=false;
end

if(nargin<4||isempty(numPoints))
    numPoints=500; 
end

if(nargin<3||isempty(gammaVal))
    gammaVal=18.8049;
end 

if(invertA)
    A=applyFunToEachMatrix(@inv,A);
end

numEllipse=size(z,2);

theta=linspace(-pi,pi,numPoints);
u=[cos(theta);
   sin(theta)];
xyPoints=zeros(2,numPoints,numEllipse);
for curEllipse=1:numEllipse
    Rxyrzz=A(1:2,1:2,curEllipse)*A(3,3,curEllipse);
    rz=A(1:2,3,curEllipse);
    rzzGamma=A(3,3,curEllipse)*gammaVal;
    zXYCur=z(1:2,curEllipse);

    for k=1:numPoints
        uCur=u(:,k);
        
        r=sqrt(rzzGamma/(uCur'*Rxyrzz*uCur-(rz'*uCur)^2));
        xyPoints(:,k,curEllipse)=zXYCur-r*uCur;
    end
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
