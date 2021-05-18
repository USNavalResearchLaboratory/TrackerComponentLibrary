%%DEMOMAGNETICCODE This file demonstrates the code for spherical harmonic
%                  synthesis for evaluating the magnetic flux. Plots of the
%                  deviation of true North from magnetic North and of the
%                  inclination at the reference epoch are computed and
%                  displayed. The Tracker Component Library (that this
%                  function is part of) must have been added to Matlab's
%                  search path for this function to work.
%
%The computation of the magnetic flux on a 256X256 grid of points is
%quite fast when the mex files have been compiled. When they have
%not been compiled, then it is can be very slow. Consequently, the file 
%CompileCLibraries should have been run prior to executing this example
%script.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('NOTE: The CompileCLibraries function should be run prior to')
disp('running this script or else the execution time will be far too slow.')

numPoints=256;
totalGridPoints=numPoints*numPoints;
lat=linspace(-90,90,numPoints)*pi/180;
lon=linspace(-180,180,numPoints)*pi/180;
[latGrid,lonGrid]=meshgrid(lat,lon);
latLonEllipse=[latGrid(:)';lonGrid(:)'];

%Convert from ellipsoidal latitudes to spherical latitudes (This assumes
%that the points are on the surface of the reference ellipsoid).
lonLatSpher=ellips2Sphere(latLonEllipse);
%Get the corresponding Cartesian points on the surface of the reference
%ellipsoid.
cartPoints=ellips2Cart([latLonEllipse;zeros(1,totalGridPoints)]);
%Append the spherical radii.
points=[sqrt(sum(cartPoints.*cartPoints,1));lonLatSpher];

%Due to precision limitations, it is possible that points having the same
%latitude will have spherical radii that differ by a value close to eps.
%This loop makes sure that points of the same latitude have the same radius
%value. That can greatly accelerate the evaluation of the spherical
%harmonic coefficients.
for curPoint=2:totalGridPoints
    if(points(3,curPoint)==points(3,curPoint-1))
        points(1,curPoint)=points(1,curPoint-1);
    end
end

%Load the coefficients for the magnetic model.
disp('Loading the coefficients for the WMM model.')
[C,S,a,c]=getWMMCoeffs();%One could also consider getIGRFCoeffs().
disp('Computing the magnetic flux vector on a grid of points.')
tic
[~,gradV]=spherHarmonicEval(C,S,points,a,c);
toc
B=-gradV;

disp('Determining the deviation of a magnetic compass from true North and the magnetic inclination angle.')
%We want to compute the deviation of the direction of B from Cartesian
%North in terms of radians East of North.
degEOfN=zeros(totalGridPoints,1);%Allocate space
I=zeros(totalGridPoints,1);%Allocate space
for curPoint=1:totalGridPoints
    %B(:,curPoint) defines the direction of magnetic North at the current
    %point.
    u=getENUAxes([latLonEllipse(:,curPoint);0]);
        
    %Given the field vector, only the components in the local tangent plane
    %matter.
    vEast=dot(B(:,curPoint),u(:,1));
    vNorth=dot(B(:,curPoint),u(:,2));
    
    %Find the angle East of North in degrees.
    degEOfN(curPoint)=atan2(vEast,vNorth)*180/pi;
    
    %Get the inclination vector using the dot product.
    I(curPoint)=(acos(dot(B(:,curPoint),u(:,3))/norm(B(:,curPoint)))-pi/2)*180/pi;
end

%If the inclination is not desired, a simple way to find the offset from
%the pointing direction of a compass from true North in degrees is directly
%using the command
%degEOfN=magHeading2Geog([latLonEllipse;zeros(1,totalGridPoints)],zeros(totalGridPoints,1))*180/pi;
%instead of the loops as given above.

%The direction of a compass in degrees East (clockwise) of North.
degEOfN=reshape(degEOfN,[numPoints,numPoints]);
%The inclination angle in degrees.
I=reshape(I,[numPoints,numPoints]);

figure(1)
clf
surf(lonGrid*180/pi,latGrid*180/pi,degEOfN,'FaceColor','interp',...
  'EdgeColor','none',...
  'FaceLighting','phong')

view(2)
axis tight

h1=xlabel('Longitude');
h2=ylabel('Latitude');
h3=title('Compass Direction Degrees East of North');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
colors=hot(256);
colormap([colors;colors(end:-1:1,:)])
colorbar

figure(2)
clf
surf(lonGrid*180/pi,latGrid*180/pi,I,'FaceColor','interp',...
  'EdgeColor','none',...
  'FaceLighting','phong')

view(2)
axis tight

h1=xlabel('Longitude');
h2=ylabel('Latitude');
h3=title('Inclination Angle in Degrees');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
colormap(jet(512))
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
