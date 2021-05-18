function latLonPts=approxProjEllipsHullLatLon(z,A,gammaVal,numPoints,invertA,a,f)
%%APPROXPROJELLIPSHULLLATLON Given the parameters of an ellipsoid in global
%       3D coordinates, use a local-tangent-plane approximation to project
%       the maximum extent of the uncertainty region onto the surface of
%       the curved Earth (approximated as a reference ellipsoid). The
%       outline of this region is returned in the latLonPts function.
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
%        a The semi-major axis of the reference ellipsoid. If this argument
%          is omitted or an empty matrix is passed, the value in
%          Constants.WGS84SemiMajorAxis is used.
%        f The flattening factor of the reference ellipsoid. If this
%          argument is omitted or an empty matrix is passed, the value in
%          Constants.WGS84Flattening is used.
%
%OUTPUTS: latLonPts A 2XnumPointsXN set of latitude-longitude points (in
%                   radians) for every ellipsoid projection.
%
%The function converts the ellipsoid into the local tangent-plane
%coordinate system at the center of the ellipse. Then, the function
%getEllipseHullPoints is called to get the projected hull on the local
%tangent plane. The points on the plane are converted back to global
%coordinates and their latitudes and longitudes are obtained. 
%
%EXAMPLE:
%This example is made with a huge uncertainty ellipse so that one can
%easily visualize its projection onto the surface of the Earth. The Earth,
%the ellipse and the projection (as a red outline) are all drawn.
% zLatLonAlt=[20.7204*(pi/180);-156.1552*(pi/180);6e3];
% zCart=ellips2Cart(zLatLonAlt);
% figure()
% clf
% hold on
% plotMapOnEllipsoid();
% A=[0.164,  -0.071,  -0.021;
%   -0.071,   0.033,   0.013;
%   -0.021,   0.013,   0.043]*1e-9;
% drawEllipse(zCart,A,[],'EdgeColor','none')
% latLonPts=approxProjEllipsHullLatLon(zCart,A);
% numPts=size(latLonPts,2);
% xyzPts=ellips2Cart([latLonPts;zeros(1,numPts)]);
% plot3(xyzPts(1,:),xyzPts(2,:),xyzPts(3,:),'-r','linewidth',4')
% view(-62,-37);
%
%December 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7||isempty(f))
    f=Constants.WGS84Flattening;
end

if(nargin<6||isempty(a))
    a=Constants.WGS84SemiMajorAxis;
end

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
latLonPts=zeros(2,numPoints,numEllipse);
for curEllipse=1:numEllipse
    plhPoint=Cart2Ellipse(z(:,curEllipse),[],a,f);
    
    %We need to establish the local tangent plane. M rotates from the
    %global coordinate system to the local tangent-plane coordinate system
    %at the center of the selected ellipse.
    M=getENUAxes(plhPoint,false,a,f)';
    
    %A holds inverse matrices we want inv(M*inv(A)*M'). This is 
    %inv(M*inv(A)*M')=inv(M')*inv(M*inv(A))=inv(M')*A*inv(M)
    %Since inv(M)=M', we then have M*A*M'.
    R=M*A(:,:,curEllipse)*M';
    %The origin of the local coordinate system is placed as
    %z(:,curEllipse).
    xyPts=getEllipseHullPoints([0;0;0],R,gammaVal,numPoints,false);
    
    %Convert the points from local tangent plane coordinates to global
    %coordinates.
    xyPts=bsxfun(@plus,M'*[xyPts;zeros(1,numPoints)],z(:,curEllipse));
    latLonAlt=Cart2Ellipse(xyPts,[],a,f);
    latLonPts(:,:,curEllipse)=latLonAlt(1:2,:);
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
