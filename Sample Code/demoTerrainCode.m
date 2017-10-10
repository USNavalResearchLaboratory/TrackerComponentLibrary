%%DEMOTERRAINCODE This file demonstrates the functions for spherical
%              harmonic synthesis of terrain data. This file requires that
%              the CompileCLibraries function has been run, because the
%              demonstration will be too slow when run in Matlab-only.
%
%This file displays plots of the shape of the Earth using the Earth2014
%model and the EGM2008 model.
%
%Note that the CompileCLibraries function should be run one before this
%this script is run.
%
%The models are run below their maximum degree and order so that each plot
%can be generated in less than 40 seconds on a 2.3GHz Intel Core i7
%processor. If more detail is desired, then the degree and order as well as
%the number of points plotted can be increased.
%
%June 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('NOTE: The CompileCLibraries function should be run prior to')
disp('running this script or else the execution time will be unusably slow.')

disp('Computing the Earth2014 Terrain Plot (with Oceans)')

%Geocentric (spherical) latitude (elevation) and longitude (azimuth) values
%in radians.
numPoints=700;
el=linspace(-90,90,numPoints)*pi/180;
az=linspace(-180,180,numPoints)*pi/180;
%The points from the meshgrid function are sorted by elevation.
[el,az]=meshgrid(el,az);
azel=[az(:)';el(:)'];

M=500;%M is the maximum degree and order of the terrain model. M<=2160.
[C,S]=getEarth2014TerrainCoeffs(M);

%Get the offsets from the ellipsoidal radius of the terrain at each
%spherical azimuth and elevation point.
terHeight=spherHarmonicEval(C,S,azel);

%The ellipsoidal height is with respect to the GRS80 ellipsoid.
a=Constants.GRS80SemiMajorAxis;
f=Constants.GRS80Flattening;
terRad=terHeight+ellipsoidalRadius(el(:),0,a,f);

%Convert the terrain points from spherical to WGS-84 ellipsoidal
%coordinates and exaggerate the ellipsoidal heights by a factor of 100.
ellipsCoord=spher2Ellipse([terRad';azel]);
ellipsCoord(3,:)=ellipsCoord(3,:)*100;

%Convert the ellipsoidal coordinates to Cartesian coordinates to plot a 3D
%Earth.
CartPoints=ellips2Cart(ellipsCoord);
XC=reshape(CartPoints(1,:),[numPoints,numPoints]);
YC=reshape(CartPoints(2,:),[numPoints,numPoints]);
ZC=reshape(CartPoints(3,:),[numPoints,numPoints]);
ZColor=reshape(ellipsCoord(3,:),[numPoints,numPoints]);

figure(1)
clf
surf(XC,YC,ZC,ZColor,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','phong')
axis square
axis tight
view([15,0])
camlight left
set(gca,'visible','off')%Do not show the axes
colormap(fliplr(hsv(256)))

disp('Computing the EGM2008 Elevation Plot (without Oceans)')

%Geocentric (spherical) latitude (elevation) and longitude (azimuth) values
%in radians.
numPoints=700;
el=linspace(-90,90,numPoints)*pi/180;
az=linspace(-180,180,numPoints)*pi/180;
%The points from the meshgrid function are sorted by elevation.
[el,az]=meshgrid(el,az);
azel=[az(:)';el(:)'];

M=600;%M is the maximum degree and order of the terrain model. M<=2160.
[C,S]=getEGM2008TerrainCoeffs(M);

%Get the orthometric heights at each spherical azimuth and elevation point.
theHeight=spherHarmonicEval(C,S,azel);

lat=spherLat2EllipsLat(el(:)');
latLon=[lat;az(:)'];

%Convert the orthometric heights to ellipsoidal heights. As the conversion
%is slow for many points and it effects will not be noticed on the scale
%given here, this line has been commented out.
%theHeight=theHeight+getEGMGeoidHeight(latLon,0,true,0);

%Convert the terrain points from spherical to WGS-84 ellipsoidal
%coordinates and exaggerate the ellipsoidal heights by a factor of 100.
ellipsCoord=[latLon;theHeight'*100];
CartPoints=ellips2Cart(ellipsCoord);
XC=reshape(CartPoints(1,:),[numPoints,numPoints]);
YC=reshape(CartPoints(2,:),[numPoints,numPoints]);
ZC=reshape(CartPoints(3,:),[numPoints,numPoints]);
ZColor=reshape(ellipsCoord(3,:),[numPoints,numPoints]);

figure(2)
clf
surf(XC,YC,ZC,ZColor,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','phong')
axis square
axis tight
view([15,0])
camlight left
set(gca,'visible','off')%Do not show the axes
colormap(fliplr(hsv(256)))

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

