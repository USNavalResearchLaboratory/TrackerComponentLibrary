function [latLonHEllipse,r,spherCent]=osculatingCoords2Ellipse(latLonHSpherPts,latLonRefEllips,aEllips,fEllips,is2D)
%%OSCULATINGCOORDS2ELLIPSE Given points in [latitude;longitude;height]
%           in the coordinates of an osculating sphere, convert it into the
%           [latitude;longitude;height] system of the reference ellipsoid.
%           defined at the point latLonRefEllips on the surface of the
%           reference ellipsoid. This is the inverse function to
%           ellips2OsculatingCoords.
%
%INPUTS: latLonHSpherPts The 3XN set of [latitude;longitude;height] points
%               to convert from osculating spherical coordinates to
%               ellipsoidal coordinates. Latitude and longitude are given
%               in radians.  If a 2XN matrix is passed, it is assumed that
%               the omitted height coordinate is zero and the height is
%               omitted from latLonHEllipse.
%   latLonRefEllips The [latitude;longitude] in radians in ellipsoidal
%               coordinates of where the osculating sphere is defined. This
%               is the point where the [latitude;longitude;height] in both
%               coordinate systems is the same (height=0).
%       aEllips The semi-major axis length of the ellipsoid. If this
%               argument is omitted or an empty matrix is passed, the value
%               in Constants.WGS84SemiMajorAxis is used.
%       fEllips The flattening factor of the ellipsoid. If this argument is
%               omitted or an empty matrix is passed, the value in
%               Constants.WGS84Flattening is used.
%          is2D Indicates whether one just cares about a 2D cut of the
%               ellipsoid. If so, then the radius used will be the radius
%               of curvature in the meridian. The default if omitted or an
%               empty matrix is passed is false.
%
%OUTPUTS: latLonHEllipse A 3XN set of points in [latitude;longitude;height]
%                    in the ellipsoidal coordinate system. If
%                    latLonHSpherPts was 2XN, then the height component of
%                    this matrix is also omitted.
%       r, spherCent The range of the opsculating sphere and the center of
%                    the sphere given in the Cartesian coordinates
%                    associated with the reference ellipse.
%
%Osculating spheres are discussed in [1]. This uses the function
%osculatingSpher4LatLon to get the parameters of the osculating sphere at
%the given reference location and then just calls specOscCoords2Ellipse.
%
%EXAMPLE:
%In this example, we plot the osculating point at a dot. We then draw the
%original ellipse. Next, we take coordinates on the surface of the
%osculating sphere, convert them to Cartesian and plot them, so that one
%can see that the osculating sphere touches the point and is not the same
%as the original ellipse, which in this case is visibly non-spherical.
% latLonRef=deg2rad([20.776659;0]);
% numPts=1000;
% latVals=linspace(-pi/2,pi/2,numPts);
% a=1;
% f=0.5;
% is2D=true;
% 
% refPt=ellips2Cart(latLonRef,a,f);
% %One half of the Earth.
% refEllipsPts1=ellips2Cart([latVals;zeros(1,numPts)],a,f);
% %The other half of the Earth.
% refEllipsPts2=ellips2Cart([latVals;pi*ones(1,numPts)],a,f);
% 
% %Convert the point of osculation to osculating coordinates and the back to
% %plot.
% latLonHOsc=ellips2OsculatingCoords([latLonRef;0],latLonRef,a,f,is2D);
% latLonHOscBack=osculatingCoords2Ellipse(latLonHOsc,latLonRef,a,f,is2D);
% refPtOscBack=ellips2Cart(latLonHOscBack,a,f);
% 
% %Take the ellipsoidal points as being in osculating points and convert then
% %to Cartesian. We get points on both sides of the cut of the osculating
% %ellipse.
% pts=[latVals;zeros(2,numPts)];
% latLonOsc1=osculatingCoords2Ellipse(pts,latLonRef,a,f,is2D);
% oscEllipsPts1=ellips2Cart(latLonOsc1,a,f);
% pts(2,:)=-pi;
% latLonOsc2=osculatingCoords2Ellipse(pts,latLonRef,a,f,is2D);
% oscEllipsPts2=ellips2Cart(latLonOsc2,a,f);
% 
% figure(1)
% clf
% hold on
% plot(refEllipsPts1(1,:),refEllipsPts1(3,:),'-k','linewidth',4)
% plot(refEllipsPts2(1,:),refEllipsPts2(3,:),'-k','linewidth',4)
% plot(oscEllipsPts1(1,:),oscEllipsPts1(3,:),'-g','linewidth',2)
% plot(oscEllipsPts2(1,:),oscEllipsPts2(3,:),'-b','linewidth',2)
% 
% scatter(refPt(1),refPt(3),500,'.r')
% scatter(refPtOscBack(1),refPtOscBack(3),100,'.b')
% axis([-2.1, 2.1, -2.1, 2.1])
%
%REFERENCES:
%[1] P. Williams and D. Last, "On Loran-C time-difference to co-ordinate
%    converters," in Proceedings of the 32nd Annual Convention & Technical
%    Symposium of the International Loran Association, Boulder, CO, 3-7
%    Nov. 2003.
%
%May 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<5||isempty(is2D))
        is2D=false;
    end
    
    if(nargin<4||isempty(fEllips))
        fEllips=Constants.WGS84Flattening;
    end

    if(nargin<3||isempty(aEllips))
        aEllips=Constants.WGS84SemiMajorAxis;
    end

    [r,spherCent]=osculatingSpher4LatLon(latLonRefEllips,aEllips,fEllips,is2D);
    latLonHEllipse=specOscCoords2Ellipse(latLonHSpherPts,r,spherCent,aEllips,fEllips);
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
