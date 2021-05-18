%%DEMOMAGNETICCOORDINATESYS This file shows how to use conversions to
%                           magnetic coordinate systems. Such coordinate
%                           systems arise when using ionospheric models.
%                           Plots of the deviation of the latitude
%                           coordinate of each system from geodetic
%                           latitude at the reference epoch are displayed.
%                           The magnetic code library (that this function
%                           is part of) must have been added to Matlab's
%                           search path for this function to work.
%
%The file CompileCLibraries should have been run prior to executing this
%example script or it will be unusably slow. The computation of apex
%coordinates and quasi-dipole coordinates can still be slow nonetheless.
%
%June 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('NOTE: The CompileCLibraries function should be run prior to')
disp('running this script or else the execution time will be far too slow.')

%A sparse grid is used so that the computations are not too slow.
numPoints=20;
totalGridPoints=numPoints*numPoints;
lat=linspace(-90,90,numPoints)*pi/180;
lon=linspace(-180,180,numPoints)*pi/180;
[latGrid,lonGrid]=meshgrid(lat,lon);
latLonEllipse=[latGrid(:)';lonGrid(:)'];

%Centered Dipole Coordinates
disp('Computing latitudes in centered-dipole coordinates')

%Convert from ellipsoidal latitudes to spherical latitudes (This assumes
%that the points are on the surface of the reference ellipsoid).
lonLatSpher=ellips2Sphere(latLonEllipse);
lonLatCD=spherITRS2SpherCD(lonLatSpher);

disp('Displaying centered-dipole latitude curves.')
figure(1)
clf
contour(lonGrid*180/pi,latGrid*180/pi,reshape(lonLatCD(2,:)*(180/pi),[numPoints,numPoints]),50)

axis tight

h1=xlabel('Longitude');
h2=ylabel('Latitude');
h3=title('CD Latitude');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
colors=jet(256);
colormap([colors;colors(end:-1:1,:)])
colorbar

%Apex corodinates
disp('Computing latitudes in magnetic apex coordinates. This can be slow.')

%Convert the positions to Cartesian coordinates; assume that the points are
%on the reference ellipsoid.
zCart=ellips2Cart([latLonEllipse;zeros(1,totalGridPoints)]);
[zApex,apexPoints,exitCode]=ITRS2MagneticApex(zCart);

disp('Displaying apex latitude curves.')
figure(2)
clf
contour(lonGrid*180/pi,latGrid*180/pi,reshape(zApex(1,:)*(180/pi),[numPoints,numPoints]),50)

axis tight

h1=xlabel('Longitude');
h2=ylabel('Latitude');
h3=title('Apex Latitude');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
colors=jet(256);
colormap([colors;colors(end:-1:1,:)])
colorbar

%Finally, we consider quasi-dipole coordinates.
disp('Computing latitudes in magnetic quasi-dipole coordinates. This can be slow.')
disp('Quasi-dipole coordinates are similar to apex coordinates.')
zCart=ellips2Cart([latLonEllipse;zeros(1,totalGridPoints)]);

zQD=ITRS2QD(zCart);

figure(3)
clf
contour(lonGrid*180/pi,latGrid*180/pi,reshape(zQD(1,:)*(180/pi),[numPoints,numPoints]),50)

axis tight

h1=xlabel('Longitude');
h2=ylabel('Latitude');
h3=title('Quasi-Dipole Latitude');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
colors=jet(256);
colormap([colors;colors(end:-1:1,:)])
colorbar

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
