function demoCurvedEarthSimSpiral()
%%DEMOCURVEDEARTHSIMSPIRAL Demonstrate the simulation of a spiraling target
%                 on a curved Earth. The target spirals along a nominal
%                 geodesic curve. For more examples of methods of
%                 determining geodesic and rhumb line (constant heading)
%                 trajectories, see the example demoCurvedEarthSimLevel.
%
%The spiraling dynamic model chosen is one where it is easy to determine
%which direction the overall flightpath will take when designing a
%trajectory.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('Spiraling flight from Mauna Loa, Hawaii to Honolulu, Hawaii at 20km initial')
disp('ellipsoidal) altitude and an average progression along the path of Mach 1')
ellipsAlt=20e3;%20km initial ellipsoidal altitude.

%The reference speed of sound for standard temperature, pressure and
%humidity.
cSoundSTP=speedOfSoundInAir();

%The Mauna Loa volcano in Hawaii.
phi=19.475*pi/180;
lambda=-155.608*pi/180;
latLonStart=[phi;lambda];
xyzStart=ellips2Cart([latLonStart;ellipsAlt]);

%Honolulu, Hawaii
phi=21.3000*pi/180;
lambda=-157.8167*pi/180;
latLonEnd=[phi;lambda];

disp('1) Computing the initial heading and distance to navigate on a geodesic curve from')
disp('Mauna Loa to Honolulu')

%This uses an approximate method of dealing with the non-zero altitude,
%which is significantly faster than the exact method.
[azStartF,distF]=indirectGeodeticProbGen(latLonStart,latLonEnd,ellipsAlt,true);

disp('3) Computing the spiraling trajectory with parameters:')
SpiralOffset=5e3;
disp(['Spiral radius:', num2str(SpiralOffset/10^3),' km'])
Nw=6;
disp(['Number of spirals along route:', num2str(Nw)])
numSteps=200;
disp(['Number of steps for the simulation:', num2str(numSteps)])

%The axes for the local ENU coordinates, which will be the initial local
%system used and thus the system in which the initial heading is specified.
uInit=getENUAxes([latLonStart;ellipsAlt]);

%When traveling the geodesic path, the total time take for the simulation.
tTotal=distF/cSoundSTP;
%The times of the state estimates.
times=linspace(0,tTotal,numSteps);
%The initial heading in the local tangent plane coordinates.
uh=[sin(azStartF);cos(azStartF);0];
vl=cSoundSTP*uh;%The linear velocity traveled in the local coordinates. 
%The turn rate to do Nw spirals over the time period allotted.
omega=2*pi*Nw/tTotal;
%The magnitude of the spiraling component of the velocity.
vsMag=omega*SpiralOffset;
%The initial spiral starts by going up (local coordinates). This should be
%orthogonal to the direction of motion.
vsInit=[0;0;1]*vsMag;
vTotalInit=vl+vsInit;
xInit=[xyzStart;vTotalInit];
aDyn=@(x,t)aSpiralSimp([x;vl;omega]);
xListGeo=RungeKCurvedAtTimes(xInit,uInit,times,aDyn);

%Convert the Cartesian locations into ellipsoidal coordinates.
latLonAlt=Cart2Ellipse(xListGeo(1:3,:));

disp('2) Plotting the trajectory and the ellipsoidal height as a function of time.')
disp('Note that the altitude does not slowly increase due to unmodeled curvature.')

figure(1)
clf
hold on
axis([-158,-158+3,19.5-1,19.5+2])
axis square
clear A

%Plot the trajectory of the aircraft with the curved-Earth model.
plot(latLonAlt(2,:)*(180/pi),latLonAlt(1,:)*(180/pi),'-r','linewidth',6);

h1=xlabel('Longitude');
h2=ylabel('Latitude');

set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
    
%Plot the altitude as a function of time.
figure(2)
clf
plot(times,latLonAlt(3,:),'-r','linewidth',6);
axis([0 times(end) 0 25e3])

h1=xlabel('Time (Seconds)');
h2=ylabel('Altitude (Meters)');

set(gca,'FontSize',18,'FontWeight','bold','FontName','Times')
set(h1,'FontSize',20,'FontWeight','bold','FontName','Times')
set(h2,'FontSize',20,'FontWeight','bold','FontName','Times')
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
