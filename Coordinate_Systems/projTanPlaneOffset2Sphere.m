function latLonProj=projTanPlaneOffset2Sphere(distEN,latLonRef,rE,useScalingApprox)
%%PROJTANPLANEOFFSET2SPHERE Given the location of where a local tangent
%       plane is positioned on a reference sphere, as well as offsets in
%       the East and North directions, project the point on the tangent
%       plane onto the reference sphere. This can be done one of two ways.
%
%INPUTS: distEN A 2XN set of N [East;North} distances in the local tangent
%               plane to travel.
%     latLonRef The [latitude;longitude] in radians of the reference point
%               defining the tangent plane. These are with respect to the
%               reference sphere.
%            rE The radius of the reference sphere. If this is omitted or
%               an empty matrix is passed, then the default of 
%               osculatingSpher4LatLon(latLonRef) is used.
% useScalingApprox If this is false, then this function just find the point
%               on the tangent plane at sets its height above the sphere to
%               0 to get the projected location. Otherwise, as described in
%               [1], an approximation is made to reduce the bias. This bias
%               reduction method is equivalent to using azimuthal
%               equidistant coordinates. The default if omitted or an empty
%               matrix is passed is false.
%
%OUTPUTS: latLonProj A 2XN set of projected [latitude;longitude] locations
%                    on the surface of the reference sphere in radians.
%
%EXAMPLE:
%In this example, a black half-circle representing points on the equator is
%drawn. A tangent plane is then shown in cyan, with the tangent point
%marked in blue. Lines are drawn between points along the tangent playing.
%The red lines show where those tangent plane points could project when
%setting useScalingApprox=false and green lines show where they would
%project when setting useScalingApprox=true.
% rE=1;
% numLinesDown=10;
% numPts=100;
% lonPts=linspace(0,pi,numPts);
% 
% latLonPts=[zeros(1,numPts);lonPts];
% CartPts=ellips2Cart(latLonPts,rE,0);
% latLonRef=[0;pi/2];
% tanPtCart=ellips2Cart(latLonRef,rE,0);
% 
% %Draw the points on the equator.
% figure(1)
% clf
% hold on
% %Draw half the equator.
% plot(CartPts(1,:),CartPts(2,:),'-k','linewidth',2)
% %The the slice of the tangent plane, so we just need to draw the +/-East
% %direction.
% uENU=getENUAxes(latLonRef,false,rE,0);
% uE=uENU(:,1);
% startPt=tanPtCart-uE;
% endPt=tanPtCart+uE;
% plot([startPt(1),endPt(1)],[startPt(2);endPt(2)],'-c','linewidth',2)
% 
% %For a number of points in that tangent plane, draw a line to where it maps
% %using each method.
% distEast=linspace(-1,1,numLinesDown);
% for k=1:numLinesDown
%     pointInPlane=tanPtCart+uE*distEast(k);
%     scatter(pointInPlane(1),pointInPlane(2),200,'.r','linewidth',2)
%     useScalingApprox=false;
%     latLonProj=projTanPlaneOffset2Sphere([distEast(k);0],latLonRef,rE,useScalingApprox);
%     xCart=ellips2Cart(latLonProj,rE,0);
%     plot([pointInPlane(1);xCart(1)],[pointInPlane(2);xCart(2)],'-r','linewidth',2)
% 
%     useScalingApprox=true;
%     latLonProj=projTanPlaneOffset2Sphere([distEast(k);0],latLonRef,rE,useScalingApprox);
%     xCart=ellips2Cart(latLonProj,rE,0);
%     plot([pointInPlane(1);xCart(1)],[pointInPlane(2);xCart(2)],'-g','linewidth',2)
% end
% %Mark the tangent point.
% scatter(tanPtCart(1),tanPtCart(2),400,'.b')
%
%REFERENCES:
%[1] David F. Crouse, "Worldwide ground target state propagation," in
%    Proceedings of the 23rd International Conference on Information
%    Fusion, Virtual Conference, 6-9 Jul. 2020.
%
%April 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(useScalingApprox))
    useScalingApprox=false;
end

if(nargin<3||isempty(rE))
    rE=osculatingSpher4LatLon(latLonRef);
end

if(useScalingApprox)
    dLocal=sqrt(sum(distEN.^2,1));
    distEN=bsxfun(@rdivide,distEN,dLocal);%Normalize.
    dTan=rE*tan(dLocal/rE);
    distEN=bsxfun(@times,dTan,distEN);%Rescale
end

u=getENUAxes(latLonRef,false,rE,0);
uEast=u(:,1);
uNorth=u(:,2);

cartLocRef=ellips2Cart(latLonRef,rE,0);
%The reference point plus the offsets.
pointsCart=bsxfun(@plus,cartLocRef,bsxfun(@times,distEN(1,:),uEast)+bsxfun(@times,distEN(2,:),uNorth));
latLonProj=Cart2Ellipse(pointsCart,[],rE,0);
latLonProj=latLonProj(1:2,:);

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
