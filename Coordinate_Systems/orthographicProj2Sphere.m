function latLonPts=orthographicProj2Sphere(xy,latLonRef,r)
%%ORTHOGRAPHICPROJ2SPHERE Given a point in a tangent plane with respect to
%           a point on the sphere, project the point onto the sphere using
%           an orthographic projection. A line from the point in the plane
%           to the sphere below the plane is normal to the plane (not the
%           sphere).
%
%INPUTS: xy A 2XN set of coordinates of the points in the tangent plane
%           defined at latLonRef that should be projected onto the
%           reference sphere.  The 3D location of the ith point in the
%           plane is xyPts(1,i)*uEast+xyPts(2,i)*uNorth, where uEast and
%           uNorth are the first two vectors in the output of
%           getENUAxes(latLonRef,false,r,0). Note that for each point
%           norm(xy(:,i))<=r to be valid.
% latLonRef The 2X1 [latitude;longitude] reference point for the
%           orthographic projection.
%         r The radius of the reference sphere. If this is omitted or an
%           empty matrix is passed, then osculatingSpher4LatLon(latLonRef)
%           is used.
%
%OUTPUTS: latLonPts The 2XN set of the points in xy projected
%                   orthographically onto the reference sphere.
%
%The conversion from Chapter 20 of [1] is implemented.
%
%EXAMPLE:
%In this example, we plot an x-z cut of a sphere. The reference point is
%taken to lie in that cut. We choose points on the tangent plane that lie
%within the cut and plot them. We then get their projections onto the
%sphere and plot the projections and draw lines between the points and
%their projections. That way, one can see that the lines hit the plane at
%right angles.
% r=1;
% numPts=1000;
% latLonRef=[pi/4;0];
% latLonRefCart=ellips2Cart(latLonRef,r,0);
% latVals=linspace(-pi/2,pi/2,numPts);
% spherPoints1=ellips2Cart([latVals;zeros(1,numPts)],r,0);
% spherPoints2=ellips2Cart([latVals;pi*ones(1,numPts)],r,0);
% 
% %Select some point in the plane and find their projections on the sphere.
% numYPts=11;
% xyPts=[zeros(1,numYPts);linspace(-r,r,numYPts)];
% latLonPts=orthographicProj2Sphere(xyPts,latLonRef,r);
% ptsOnSpher=ellips2Cart(latLonPts,r,0);
% 
% %Take those points on the plane and convert them to unprojected global
% %coordinates.
% uENU=getENUAxes(latLonRef,false,r,0);
% xyPts3D=latLonRefCart+bsxfun(@times,xyPts(1,:),uENU(:,1))+bsxfun(@times,xyPts(2,:),uENU(:,2));
% 
% figure(1)
% clf
% hold on
% plot(spherPoints1(1,:),spherPoints1(3,:),'-k','linewidth',2)
% plot(spherPoints2(1,:),spherPoints2(3,:),'-k','linewidth',2)
% %Draw a line representing the tangent plane by connecting the first and
% %last of the points in the plane.
% plot([xyPts3D(1,1),xyPts3D(1,end)],[xyPts3D(3,1),xyPts3D(3,end)],'-b','linewidth',2)
% 
% %Draw lines between the projections points and the points on the sphere.
% for k=1:numYPts
%      plot([xyPts3D(1,k),ptsOnSpher(1,k)],[xyPts3D(3,k),ptsOnSpher(3,k)],'-c')
% end
% scatter(xyPts3D(1,:),xyPts3D(3,:),200,'.r')
% scatter(ptsOnSpher(1,:),ptsOnSpher(3,:),200,'.r')
% axis square
% axis equal
%
%REFERENCES:
%[1] J. P. Snyder, "Map projections- a working manual," U.S. Geological
%    Survey, Tech. Rep. 1395, 1987.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(r))
    r=osculatingSpher4LatLon(latLonRef); 
end

%Sine and cosine of the reference latitude.
sinPhi1=sin(latLonRef(1));
cosPhi1=cos(latLonRef(1));
%Reference longitude.
lambda0=latLonRef(2);

x=xy(1,:);
%We are 
y=xy(2,:);

%Equation 20-18.
rho=sqrt(x.^2+y.^2);
%Equation 20-19, removing the asin. Also, the min is to deal with points
%that due to finite precision errors might be just a little bit larger in
%magnitude than r even though they should really be r.
sinC=min(1,rho./r);
%Taking the cosine of Equation 20-19 and choosing the positive solution. To
cosC=sqrt(1-sinC.^2);

%Equation 20-14 and 20-15.
latLonPts=[asin(cosC*sinPhi1+y.*sinC.*cosPhi1./rho);
           lambda0+atan2(x.*sinC,rho.*cosPhi1.*cosC-y.*sinPhi1.*sinC)];
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
