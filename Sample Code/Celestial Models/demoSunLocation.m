function demoSunLocation()
%%DEMOSUNLOCATION Demonstrate how the location of the sun can be found. The
%                 apparent locations of the Sun with and without the
%                 effects of atmospheric refraction (using a very simple
%                 refraction model) are displayed for an observer on the
%                 reference ellipsoid at Hilo, Hawaii on 1 June 2013.
%                 The apparent outline of the Sun is shown. The accuracy of
%                 the image varies largely based on atmospheric turbulence
%                 and the difficulty in accurately modeling atmospheric
%                 refraction at low elevation angles.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('Simulating the location of the Sun with and without refraction for an')
disp('Observer in Hilo, Hawaii on 1 June 2013 at 15:42 UTC.')
%The Latitude and Longitude of Hilo, Hawaii.
lat=19.7056*(pi/180);
lon=-155.0858*(pi/180);
%An observer at sea level.
obsLoc=ellips2Cart([lat;lon;0]);

%On June 1st 2013, the apparent sunrise in Honolulu was at about 15:41 UTC.
%Thus, we choose 15:42UTC, to be about halfway through the sunrise.
[JulUTC1,JulUTC2]=Cal2UTC(2013,6,1,15,42,0);

disp(' ')%Add a line break.
disp('1) Determining the refraction-free outline of the Sun')

%To display an image of the outline of the Sun, the radius of the Sun is
%needed. It is (in meters).
sunRad=696000e3;

%Find the Cartesian location of the Sun in WGS-84 ECEF coordinates,
%correcting for aberration and light-time, but not for atmospheric
%refraction.
[~,rSunITRS]=solarBodyVec(JulUTC1,JulUTC2,'UTC',11,[obsLoc;0;0;0],'ITRS');
rSunITRS=rSunITRS(1:3);%Only keep the position.

%The location of the Sun with respect to the observer.
rSunLoc=rSunITRS-obsLoc;

%A rotation matrix that rotates the z-axis to the direction of
%the Sun.
R=rotAxis2Vec(rSunLoc/norm(rSunLoc));

%Rotate a set of unit vectors such that x-and y are in the
%plane orthogonal to the vector going from the observer to the
%Sun.
uBasis=R*eye(3);

%We will find the outer points of the Sun that are seen by
%the observer.
numPoints=50;
polAng=linspace(0,2*pi,numPoints);
%The locations of the outline points of the Sun in ITRS (A circle).
outlineLocs=bsxfun(@plus,rSunLoc,sunRad*bsxfun(@times,cos(polAng),uBasis(:,1))+sunRad*bsxfun(@times,sin(polAng),uBasis(:,2)));

%This will be converted to local ENU coordinates for plotting. First, we
%need to determine the relationship between ECEF and ENU coordinates.
%Get the ENU axes.
uENU=getENUAxes([lat;lon;0]);

%Find the rotation matrix to go from ECEF to ENU.
ECEF2ENURotMat=findTransParam(eye(3),uENU);

%Rotate into the local ENU coordinate system
for curPoint=1:numPoints
    outlineLocs(:,curPoint)=ECEF2ENURotMat*outlineLocs(:,curPoint);
end
%Get the azimuth and elevation of the Sun in spherical coordinates with
%respect to local ENU coordinates.
outlineSphere=Cart2Sphere(outlineLocs);
azSun=outlineSphere(2,:);
elSun=outlineSphere(3,:);

disp(' ')%Add a line break.
disp('2) Determining the refraction-corrupted outline of the Sun, assuming yellow light.')
disp('In simple refraction models, refraction only changes the elevation of the')
disp('Sun in the sky. This step is slow, because it uses ray tracing.')
%Add in the effects of standard refraction for yellow light using the
%simple ray tracing algorithm. Say 20 degrees Centigrade, 293.15K, 70%
%humidity and standard pressure.
zenithDist=pi/2-elSun;
apparentZenithDist=addAstroRefrac(0,[lat;lon;0],zenithDist,0.7,Constants.standardAtmosphericPressure,293.15);
apparentElSun=pi/2-apparentZenithDist;

disp(' ')%Add a line break.
disp('3) Plotting the outlines of the Sun. The refraction-corrupted image is')
disp('the higher one.')

%For determining good bounds for the plot, the azimuth and elevation of the
%Sun in the local ENU coordinate system need to be found. 
%First, get a unit vector

uSun=rSunITRS/norm(rSunITRS);

%The location of the center of the sun in ENU coordinates.
uSunENU=ECEF2ENURotMat*uSun;

%Spherical coordinate location of the Sun in ENU.
sunSphere=Cart2Sphere(uSunENU);
sunAz=sunSphere(2);
sunEl=sunSphere(3);

figure(1)
clf
hold on
%Plot the refraction-free Sun.
plot(90-azSun*(180/pi),elSun*(180/pi),'linewidth',2)

%Plot the refraction-corrupted Sun
plot(90-azSun*(180/pi),apparentElSun*(180/pi),'-r','linewidth',4)

%Plot the apparent horizon within +/-2 degrees of the Sun location.
polAng=linspace((pi/2-sunAz)-1*(pi/180),(pi/2-sunAz)+1*(pi/180),numPoints);
plot(polAng*(180/pi),zeros(1,numPoints),'-k','linewidth',2)

axis([90-sunAz*(180/pi)-1,90-sunAz*(180/pi)+1,sunEl*(180/pi)-1,sunEl*(180/pi)+1])
axis square
h1=xlabel('Degrees East of North');
h2=ylabel('Degrees Elevation');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')

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
