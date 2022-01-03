function h=drawProjEllipsHull(z,A,gammaVal,AIsInv,varargin)
%%DRAWPROJELLIPSHULL Consider an ellipsoid in 3D (x,y,z) coordinates. For
%       every fixed z coordinate, if the x-y plane at the z coordinate
%       intersects the ellipsoid, it cuts an ellipse. For example, the
%       function projEllipse2ZPlane gets the parameters of the ellipse of
%       intersection. This function plots the lines marking the outer
%       limits of all of those ellipses in the x-y plane (the hull of the
%       ellipses). This can be used for plotting the maximum uncertainty in
%       x-y (say, on the ground) irrespective of the z coordinate (maybe
%       the altitude).
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
%  invertA If this is true, then A is inverted before use. The default if
%          omitted or an empty matrix is passed is false.
% varargin Sets of values that should be passed to the plot function to
%          format the ellipses or that will be passed to the plot function
%          to format the drawn line. For example, 'linewidth', 2 will make
%          the line thicker.
%
%OUTPUTS: h An NX1 cell array containing the plot objects for each of the
%           ellipses. This can be useful if, for example, one wishes to
%           change the transparency of object i to 50%, one can use the
%           alpha(h{i},0.5) command.
%
%This function just calls getEllipseHullPoints to get a number of points on
%the hull and then plots the result.
%
%EXAMPLE:
%This creates a 3D ellipsoid with a random rotation, plots the ellipsoid in
%3D and then plots the 2D hull projection around it. Run it a few times and
%one can see how the hull projection marks the maximum extent of the
%ellipsoid in x-y across all z values.
% eigVals=[20;8;1];
% rotMat=randRotMat(3);
% A=rotMat*diag(eigVals)*rotMat';
% z=[10;0;0];
% figure(1)
% clf
% hold on
% drawEllipse(z,A)
% drawProjEllipsHull(z,A,[],[],'linewidth',2)
% view(-10,11)
%
%December 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(AIsInv))
    AIsInv=true;
end

if(nargin<3||isempty(gammaVal))
    gammaVal=18.8049;
end

if(AIsInv==false)
    A=applyFunToEachMatrix(@inv,A);
end

%Save the value of hold on the plot so that it can be reset to its previous
%state after plotting all of the ellipses.
holdVal=ishold();

numEllips=size(z,2);
numPoints=500;
xyPoints=getEllipseHullPoints(z,A,gammaVal,numPoints,false);
for curEllipse=1:numEllips    
    h=plot(xyPoints(1,:),xyPoints(2,:),varargin{:});
end

%Restore the hold value to its original setting.
if(~holdVal)
    hold off
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
