%%DEMOGRAVCODEMOONMODEL This file demonstrates how to load the spherical
%                       harmonic coefficients associated with the Moon. It
%                       then makes a plot of the gravitational anomaly and
%                       the deflection of the vertical compared to an
%                       ellipsoidal model for points 2km above the surface.
%
%Note that the CompileCLibraries function only needs to be run the first
%time this script is run. It can be commented out on subsequent runs.
%
%Note that the CompileCLibraries function should be run one before this
%this script is run.
%
%The models are run below their maximum degree and order so that the plots
%can be rapidly generated.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

display('NOTE: The CompileCLibraries function should be run prior to')
display('running this script or else the execution time will be unusably slow.')

%%%%%PLOT THE TIDE-FREE GRAVITY ANOMALY AND VERTICAL DEFLECTION AT 2KM ELEVATION%%%%%
display('Computing the Tide-Free Gravity Anomaly and Vertical Deflection Plots')
numPoints=500;
%Geodetic (ellipsoidal) latitude and longitude in radians.
lat=linspace(-90,90,numPoints)*pi/180;
lon=linspace(-180,180,numPoints)*pi/180;
[latGrid,lonGrid]=meshgrid(lat,lon);
numGridPoints=numPoints^2;
%Note that the latLon values were entered into the meshgrid function in
%such a way that the points in latLon are sorted by latitude. This allows
%the spherical harmonic synthesis to be greatly accelerated.
latLon=[latGrid(:)';lonGrid(:)'];

%The maximum degree and order of the coefficients for the geopotential
%model. M<=900 for the GL0900C model.
M=360;
%Use coefficients for the tide-free model rather than for the zero-tide
%model.
isTideFree=true;
%Get the spherical harmonic coefficients for the lunar gravitational model.
[C,S,a,c]=getMoonGravCoeffs(M,isTideFree);

%Compute the normal gravity 2km (2e3m) above the reference ellipsoid at all
%of the given ellipsoidal latitudes and longitudes. The JPL parameters for
%a cooridnate system on the Moon are used.
aMoon=Constants.JPLMoonSemiMajorAxis;
fMoon=Constants.JPLMoonFlattening;

cartPoints=ellips2Cart([latGrid(:)';lonGrid(:)';9e3*ones(1,numGridPoints)],aMoon,fMoon);
%To have an ellipsoidal gravity model, a value of GM is needed. The value
%from the gravitational model is used.
GMMoon=c;
omegaMoon=Constants.JPLMoonRotationRate;
[~,gNormal]=ellipsParam2Grav(cartPoints,omegaMoon,aMoon,fMoon,GMMoon);

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
%acceleration of the Moon.
g=gradV+omegaMoon^2*[cartPoints(1,:);cartPoints(2,:);zeros(1,numGridPoints)];

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

h0=title('Gravity Anomaly at 2km Altitude');
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

h0=title('Deflection of the Vertical at 2km Altitude');
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
