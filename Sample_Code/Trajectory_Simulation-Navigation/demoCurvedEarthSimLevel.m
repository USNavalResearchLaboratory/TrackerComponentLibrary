function demoCurvedEarthSimLevel()
%%DEMOCURVEDEARTHSIMLEVEL  Demonstrate the simulation of a target on a
%                     curved-Earth using a simple level-flight flat-Earth
%                     target dynamic model with the WGS-84 ellipsoidal
%                     Earth model. Geodesic and rhumb-line trajectories are
%                     determined. The accuracy of a simple approximation
%                     for handling non-zero ellipsoidal altitudes is
%                     compared to slow, exact solutions. The trajectories
%                     are then plotted over 1000 points, where it is also
%                     demonstrated how the geodesic coordinate system
%                     evolution can be determined using explicit or
%                     numerically-determined derivatives.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Simulate level-flight from Hilo Hawaii to F�ssen, Germany at a constant
%altitude on an ellipsoidal Earth, traveling at a constant speed.

disp('Level Flight from Hilo, Hawaii to F�ssen, Germany at 8km (ellipsoidal) Altitude and Mach 1')
ellipsAlt=8e3;%8km ellipsoidal altitude.

%The target should be traveling at Mach 1-> the speed of sound. The speed
%of sound in air depends on the temperature, humidity, pressure, and
%composition of the air. If one were to use
cSoundSTP=speedOfSoundInAir();
%Then one would obtain a reference speed of sound for standard temperature,
%pressure and humidity.

disp(['The speed of sound at standard temperature, pressure (STP) and humidity is ',num2str(cSoundSTP),'m/s,'])
disp('which will be used to define Mach 1.')
disp(' ')%All a line-break to the output.

%The latitude and longitude of the Mauna Kea Observatory in Hilo, Hawaii.
phi=19.823*pi/180;%North latitude in radians.
lambda=-155.470*pi/180;%East longitude in radians.
latLonStart=[phi;lambda];

%The approximate latitude and longitude of Neuschwanstein Castle in F�ssen,
%Germany.
phi=47.5575*pi/180;%North latitude
lambda=10.7500*pi/180;%East longitude
latLonEnd=[phi;lambda];

disp('1) Computing the initial heading and distance to navigate on a geodesic curve from Hilo to F�ssen')
%The initial heading (in local East-North-Up coordinate) and distance
%traveled to get from Hilo, Hawaii to F�ssen, Germany is the solution to
%the indirect geodetic problem. When the altitude is not zero, traditional
%indirect geodetic problem solution methods will report a distance that is
%too small, since they assume altitude==0. Here, we will compare the
%results when assuming a zero altitude, when using an approximation with a
%non-zero altitude, and when using a more exact formula with a non-zero
%altitude. Using a non-zero altitude is significantly slower.

tic;%True means use the approximation for non-zero altitude flight.
[azStart0,dist0,azEnd0]=indirectGeodeticProb(latLonStart,latLonEnd);
time0=toc;
disp(['Assuming a zero altitude, the computation took ', num2str(time0), ' seconds.'])

tic;%True means use the approximation for non-zero altitude flight.
[azStartF,distF]=indirectGeodeticProbGen(latLonStart,latLonEnd,ellipsAlt,true);
time1=toc;

disp(['Using an approximation to handle the non-zero altitude took ', num2str(time1), ' seconds.'])

tic;%False means use the exact iteration for non-zero altitude flight.
[azStartE,distE]=indirectGeodeticProbGen(latLonStart,latLonEnd,ellipsAlt,false);
time2=toc;
disp(['Using an exact algorithm to handle the non-zero altitude took ', num2str(time2), ' seconds.'])

disp(['The difference in distance computed between the approximation and the exact algorithms'])
disp(['is only ', num2str(abs(distE-distF)), ' meters.'])

disp(['The difference in distance computed between zero and exact non-zero altitudes is ',num2str(abs(distE-dist0)), ' meters'])
disp(['out of a total transit distance of ', num2str(distE), ' meters.'])

disp(' ')%Insert line break
disp('2) Computing the heading and distance to navigate on a rhumb line (constant-heading trajectory)')
disp('from Hilo to F�ssen')
tic;
[azimuthRhumb0, distRhumb0]=indirectRhumbProblem(latLonStart,latLonEnd);
time0=toc;
disp(['Using a zero altitude, the computation took ', num2str(time0), ' seconds.'])

tic;
[azimuthRhumbA, distRhumbA]=indirectRhumbProblem(latLonStart,latLonEnd,ellipsAlt,true);
time1=toc;
disp(['Using an approximation to handle the non-zero altitude took ', num2str(time1), ' seconds.'])

tic;
[azimuthRhumbE, distRhumbE]=indirectRhumbProblem(latLonStart,latLonEnd,ellipsAlt,false);
time2=toc;
disp(['Using an exact algorithm to handle the non-zero altitude took ', num2str(time2), ' seconds.'])

disp(['The difference in distance computed between the approximation and the exact algorithms is ', num2str(abs(distRhumbE-distRhumbA)), ' meters.'])

disp(['The difference in distance computed between zero and exact non-zero altitudes is ',num2str(abs(distRhumbE-distRhumb0)), ' meters'])
disp(['out of a total transit distance of ', num2str(distRhumbE), ' meters.'])

disp(' ')%Insert line break
disp('3) Computing the trajectory at 1000 points of a non-maneuvering target flying the geodesic path using')
disp('analytical derivatives.')
%The continuous-time drift function for a 3D, flat-Earth dynamic model is
aFlatEarth3D=@(x,t)aPoly(x,3);
%The initial local coordinate system is defined as East-North-Up,
%corresponding to the x, y and z axes. This lets us easily translate the
%initial heading given by the indirectGeodeticProb function into components
%for the initial heading. The basis vectors at the starting point are
uInit=getENUAxes(latLonStart);

%Obtain the initial Cartesian location of the target.
rGlobalCart=ellips2Cart([latLonStart;ellipsAlt]);
%The initial velocity in LOCAL coordinates of the target is Mach 1 times
%the direction obtained by solving the indirect geodesic problem. The
%direction is found from trigonometry using the heading (in radians East of
%North).
vLocalInit=cSoundSTP*[sin(azStartE);cos(azStartE);0];

%The initial target state with GLOBAL position and LOCAL velocity is thus
xInit=[rGlobalCart;vLocalInit];

%When traveling in a straight line on an ellipsoidal Earth, the local basis
%vectors evolve according to
uDyn=@(u,x,t)uDotEllipsoid(u,x);

%The target is going to travel a distance of distE meters and it is
%traveling a constant cSoundSTP meters per second. Thus, the total travel
%time of the target is
tTotal=distE/cSoundSTP;

%Determine the target location at 1000 times along the trajectory. Since we
%are not setting a minimum Runge-Kutta step size, this will also be where
%all of the Runge-Kutta integration steps are placed.
times=linspace(0,tTotal,1000);

%Perform Runge-Kutta integration to determine the location of the target.
%The function RungeKCurvedAtTimes maps the flat-Earth model to a curved
%Earth.
tic;
xList=RungeKCurvedAtTimes(xInit,uInit,times,aFlatEarth3D,uDyn);
timeGeo=toc;

disp(['Analytic trajectory determination took ', num2str(timeGeo),' seconds.'])

disp('Computing the trajectory at 1000 points of a non-maneuvering target flying the geodesic path using')
disp('numerical differentiation.')
%For travel on/ over a general surface, where determining the
%time-derivatives of the local basis vectors can be difficult, one can use
%a numerical differentiation solution as long as one can express the local
%vertical anywhere on/ above the surface (curvature of the plumb line is
%not tolerated). Here, the numerical solution algorithm is demonstrated
%just using the ellipsoidal Earth model.

downVecFunc=@(r,t)(-getENUAxes(Cart2Ellipse(r),true));
uDer=@(u,x,t)uDotNumeric(u,x,t,downVecFunc);
tic;
xListNumeric=RungeKCurvedAtTimes(xInit,uInit,times,aFlatEarth3D,uDer);
timeGeoNum=toc;
disp(['Numerical trajectory determination took ', num2str(timeGeoNum),' seconds.'])

disp('Computing the trajectory at 1000 points of a non-maneuvering target flying the rhumb path.')

%The initial heading for the Rhumb-line path is different than that of the
%geodesic curve, because it does not change with time.
vLocalInit=cSoundSTP*[sin(azimuthRhumbE);cos(azimuthRhumbE);0];

%The initial target state with GLOBAL position and LOCAL velocity is thus
xInit=[rGlobalCart;vLocalInit];

%When traveling along a rhumb-line, the local coordinate system is known at
%all times: it is the local ENU coordinate system. The local basis vectors
%are deterministially known at all places on the Earth.
uDet=@(x,t)getENUAxes(Cart2Ellipse(x(1:3)));

%The target is going to travel a distance of distRhumb meters and it is
%traveling a constant cSoundSTP meters per second. Thus, the total travel
%time of the target is
tTotal=distRhumbE/cSoundSTP;

%Determine the target location at 1000 times along the trajectory. Since we
%are not setting a minimum Runge-Kutta step size, this will also be where
%all of the Runge-Kutta integration steps are placed.
times=linspace(0,tTotal,1000);

%Perform Runge-Kutta integration to determine the location of the target.
%The function RungeKCurvedAtTimes maps the flat-Earth model to a curved
%Earth. The lack of an initial uInit tells the function that uDet provides
%deterministic basis vectors and not basis vector derivatives.
tic
xListRhumb=RungeKCurvedAtTimes(xInit,[],times,aFlatEarth3D,uDet);
timeRhumb=toc;

disp(['Rhumb trajectory determination took ', num2str(timeRhumb),' seconds.'])

disp(' ')%Insert line break
disp('4) Displaying The Geodesic and Rhumb Trajectories. Red=Geodesic; Green=Rhumb.')

%Display a map of the Earth
figure(1)
clf
plotMapOnEllipsoid();
view(22,59);%Specified the viewpoint of the camera.
hold on

%Plot the simulated trajectory on the Earth. 
plot3(xList(1,:),xList(2,:),xList(3,:),'-r','linewidth',9)

%Plot the simulated rhumb trajectory on the Earth. 
plot3(xListRhumb(1,:),xListRhumb(2,:),xListRhumb(3,:),'-g','linewidth',9)

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
