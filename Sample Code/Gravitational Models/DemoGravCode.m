%%DEMOGRAVCODE This file demonstrates the functions for spherical harmonic
%              synthesis of gravitational data. This file requires that the
%              CompileCLibraries function has been run, because the
%              demonstration will be too slow when run in Matlab-only. Use
%              the function DemoGravCodeNoCompile for an example that can
%              be run without compiling the libraries. 
%
%This file displays a plot of the geoid using the EGM2008 model,
%exaggerating the geoid heights by a factor of 10^4. It then plots the
%tide-free gravity anomaly and vertical deflection at 9km altitude
%(ignoring the atmosphere) using the EGM2008 model.
%
%Note that the CompileCLibraries function should be run one before this
%this script is run.
%
%The models are run below their maximum degree and order so that each plot
%can be generated in less than 40 seconds on a 2.3GHz Intel Core i7
%processor. If more detail is desired, then the degree and order as well as
%the number of points plotted can be increased.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('NOTE: The CompileCLibraries function should be run prior to')
disp('running this script or else the execution time will be unusably slow.')

%%%%%PLOT THE GEOID HEIGHT%%%%%
disp('Computing the Geoid Height Plot')
numPoints=700;
%Geodetic (ellipsoidal) latitude and longitude values in radians.
lat=linspace(-90,90,numPoints)*pi/180;
lon=linspace(-180,180,numPoints)*pi/180;
[latGrid,lonGrid]=meshgrid(lat,lon);
%Note that the latLon values were entered into the meshgrid function in
%such a way that the points in latLon are sorted by latitude. This allows
%the spherical harmonic synthesis to be greatly accelerated.
latLon=[latGrid(:)';lonGrid(:)'];

%The function getGeoidHeight will automatically load the
%coefficients for the gravitational model to the maximum degree and order
%if they are not provided. However, if multiple calls are made to the
%function, it is faster if the coefficients are pre-loaded. In this case,
%the coefficients are being preloaded so that a lower degree and order can
%be used than the maximum, which is the default.
%The maximum degree and order of the coefficients for the disturbing model
%to use. M<=2190 for the EGM2008 model.
M=360;
%The model will be used with the same approximations used by the NGA so
%that the results are consistent with theirs. See the documentation for the =
%function getEGMGeoidHeight for more details.
useNGAApprox=true;
%Get the coefficients for the disturbing potential.
[C,S]=getEGMWGS84TCoeffs(M,useNGAApprox);
coeffData.C=C;
coeffData.S=S;

MZeta=360;%The maximum degree and order of the correction term model.
%MZeta<=2160 in the EGM2008 model.
[C,S]=getEGMZeta2NCoeffs(MZeta);
coeffData.CZeta=C;
coeffData.SZeta=S;

modelType=0;%Use the EGM2008 model.
%Use a tide-free geoid model.
tideSys=0;
geoidHeight=getEGMGeoidHeight(latLon,tideSys,useNGAApprox,modelType,coeffData);

%Transform the latitude, longitude, and ellipsoidal height coordinates
%into Cartesian coordinates and plot the geoid in 3D with the height
%component exaggerated by a factor of 10^4
cartPoints=ellips2Cart([latGrid(:)';lonGrid(:)';1e4*geoidHeight(:)']);
XC=reshape(cartPoints(1,:),[numPoints,numPoints]);
YC=reshape(cartPoints(2,:),[numPoints,numPoints]);
ZC=reshape(cartPoints(3,:),[numPoints,numPoints]);
ZColor=reshape(geoidHeight,[numPoints,numPoints]);

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
colormap(jet(256))

%%%%%PLOT THE TIDE-FREE GRAVITY ANOMALY AND VERTICAL DEFLECTION AT 9KM ELEVATION%%%%%
disp('Computing the Tide-Free Gravity Anomaly and Vertical Deflection Plots')
%The gravity anomaly is the difference between the magnitude of the "true"
%acceleration due to gravity (in this case using the EGM2008 model) and
%acceleration due to gravity from a model, in this case the WGS84
%ellipsoidal Earth model. The 9km elevation above the reference ellipsoid
%was chosen as all points are above the terain. The effects of being in the
%atmosphere are ignored.

numPoints=300;
%Geodetic (ellipsoidal) latitude and longitude values in radians.
lat=linspace(-90,90,numPoints)*pi/180;
lon=linspace(-180,180,numPoints)*pi/180;
[latGrid,lonGrid]=meshgrid(lat,lon);
numGridPoints=numPoints^2;
%Note that the latLon values were entered into the meshgrid function in
%such a way that the points in latLon are sorted by latitude. This allows
%the spherical harmonic synthesis to be greatly accelerated.
latLon=[latGrid(:)';lonGrid(:)'];

%The maximum degree and order of the coefficients for the geopotential
%model. M<=2190 for the EGM2008 model.
M=360;
%Use coefficients for the tide-free model rather than for the zero-tide
%model.
isTideFree=true;
%Get the spherical harmonic coefficients for the gravitational model.
[C,S,a,c]=getEGMGravCoeffs(M,isTideFree);

%Compute the normal gravity 9km (9e3m) above the reference ellipsoid at all
%of the given ellipsoidal latitudes and longitudes.
cartPoints=ellips2Cart([latGrid(:)';lonGrid(:)';9e3*ones(1,numGridPoints)]);
[~,gNormal]=ellipsParam2Grav(cartPoints);

%Compute the acceleration due to gravity using the spherical harmonic
%coefficients. For this purpose, the points need to be converted to
%spherical coordinates.
r=sqrt(sum(cartPoints.*cartPoints,1));
spherPoints=[r;ellips2Sphere(latLon)];
%Due to precision limitations, it is possible that points having the same
%latitude will have spherical radii that differ by a value close to eps.
%This loop makes sure that points of the same latitude have the same radius
%value. That can greatly accelerate the evaluation of the spherical
%harmonic coefficients when the points are evaluated on a grid and have
%been sorted by latitude.
for curPoint=2:numGridPoints
    if(spherPoints(3,curPoint)==spherPoints(3,curPoint-1))
        spherPoints(1,curPoint)=spherPoints(1,curPoint-1);
    end
end

[~,gradV]=spherHarmonicEval(C,S,spherPoints,a,c);

%The gravity is the acceleration due to gravity plus the centripetal
%acceleration of the Earth.
omega=Constants.WGS84EarthRotationRate;
%We approximate the rotation axis of the Earth as the z-axis. Earth
%orientation parameters could be used to get a more accurate approximation.
g=gradV+omega^2*[cartPoints(1,:);cartPoints(2,:);zeros(1,numGridPoints)];

%The gravity anomaly is the difference in magnitude between the gravity and
%normal gravity.
gravAnom=reshape(sqrt(sum(g.*g,1))-sqrt(sum(gNormal.*gNormal,1)),[numPoints,numPoints]);

%The vertical deflection is the difference in angle between the gravity and
%normal gravity. The dot product is used to extract the angle between the
%vectors.
vertDeflect=reshape(acos(sum(gNormal.*g,1)./(sqrt(sum(gNormal.*gNormal,1)).*sqrt(sum(g.*g,1)))),[numPoints,numPoints]);

lSize=14;
fSize=14;

figure(2)
clf
%The factor of 1e-5 converts the anomalies into milliGals.
surf(lonGrid*180/pi,latGrid*180/pi,gravAnom/1e-5,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','phong')
view(2)
axis tight

h0=title('Gravity Anomaly at 9km Altitude');
h1=xlabel('Longitude');
h2=ylabel('Latitude');
set(gca,'FontSize',lSize,'FontWeight','bold','FontName','Times')
set(h0,'FontSize',fSize,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',fSize,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',fSize,'FontWeight','bold','FontName','Times')
colormap(gray(256))
set(get(colorbar('peer',gca),'ylabel'),'string','mGal','FontSize',lSize,'FontWeight','bold','FontName','Times')

figure(3)
clf
%The factor of 180/pi converts the radians into degrees. The 60*60
%converts that to arcseconds.
surf(lonGrid*180/pi,latGrid*180/pi,vertDeflect*(180/pi)*60*60,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','phong')
view(2)
axis tight

h0=title('Deflection of the Vertical at 9km Altitude');
h1=xlabel('Longitude');
h2=ylabel('Latitude');
set(gca,'FontSize',lSize,'FontWeight','bold','FontName','Times')
set(h0,'FontSize',fSize,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',fSize,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',fSize,'FontWeight','bold','FontName','Times')
colormap(jet(256))
set(get(colorbar('peer',gca),'ylabel'),'string','Arcseconds','FontSize',lSize,'FontWeight','bold','FontName','Times')

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
