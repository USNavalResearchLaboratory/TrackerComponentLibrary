function demoStarTrails()
%%DEMOSTARTRAILS This demonstrates how to determine the apparent locations
%                of stars in the sky at a particular location and time
%                using the Hipparcos catalog and a low-fidelity atmospheric
%                refraction model. The stars are plotted over a number of
%                times on the same plot such that their motion in the sky
%                becomes apparent. The planets, Sun and Moon are not
%                plotted.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The location of the observer on the surface of the Earth in ellipsoidal
%coordinates.
obsLoc=[-pi/8;pi-pi/8+pi/8;0];

%Properties of the imaging sensor:
pixX=3264;%Pixels wide
pixY=2448;%Pixels high.
pixWidth=1.5e-6;%Width of a pixel in meters
f=4.12e-3;%Focal length, meters.
%The size of the sensor.
sensX=pixX*pixWidth;
sensY=pixY*pixWidth;
%The angular width of the field of view when viewing something at infinity.
radAz=2*atan(sensX/(2*f));%Width in azimuth (radians)
radEl=2*atan(sensY/(2*f));%Width in elevation (radians)

%The half-azimuth and half-elevation of the observer's field of view from
%the center of the field.
locFOVHalfWidth=[radAz/2;radEl/2];

%It is assumed that the observer's true orientation is a rotation 60
%degrees about the local ellipsoidal ENU x-axis-up followed by a rotation
%of 15 degrees about the rotated z-axis. So, we will now build the
%ENU2Local rotation matrix.
R1=axisAng2RotMat([1;0;0],60*pi/180);
R2=axisAng2RotMat([0;0;1],15*pi/180);
ENU2LocRotMat=R2*R1;

disp('Loading the star catalog.')
%Load stars with a magnitude of at most 6 (the visible region). The stars
%are formatted such that the data file can be directly passed to the
%changeEpoch function, which allows the data to be used at different times.
%The catalog is provided at epoch 1991.25.
dataFile=getHipparcosCat([],6,false,true);

%Change the epoch from 1995.25 to J2000.0
dataFile=changeEpoch(dataFile);

%Set the initial date of display. Minutes will be advanced to demonstrate
%the motion of the stars.
year=2014;
month=1;
day=1;
hour=0;
second=0;

disp('Plotting the stars at different times in the local coordinate')
disp('system and field of view of the observer')
figure(1)
hold on
for minute=0:5:20
    [Jul1,Jul2]=Cal2UTC(year,month,day,hour,minute,second);
    %This returns the apparent star locations in the local ENU coordinate
    %system, where zSpherENU is the location in spherical azimuth and
    %elevation in ENU coordinate and uObs is unit vectors in local ENU
    %coordinates.
    [zSpherENU,uStarsENU]=starCat2Obs(dataFile,Jul1,Jul2,obsLoc);

    %We have to determine which stars are in the observer's field of view.
    %First, we need to rotate the apparent locations into the observer's
    %local coordinate system
    uStarLocal=ENU2LocRotMat*uStarsENU;
    
    %Next, obtain the azimuth and elevation of the stars in the local
    %coordinate system.
    zStarLocal=Cart2Sphere(uStarLocal);
    %Discard the useless range component.
    zStarLocal=zStarLocal(2:3,:);
    
    %The thing in the middle of the field of view is at local azimuth and
    %elevation (0,0). We need to figure out which stars and celestial
    %bodies are in the true field of view and get rid of the others.
    sel=zStarLocal(1,:)<locFOVHalfWidth(1);
    sel=sel&(zStarLocal(1,:)>-locFOVHalfWidth(1));
    sel=sel&(zStarLocal(2,:)<locFOVHalfWidth(2));
    sel=sel&(zStarLocal(2,:)>-locFOVHalfWidth(2));

    uStarLocal=uStarLocal(:,sel);
    zStarLocal=zStarLocal(:,sel);

    %Plot the stars in terms of azimuth and elevation in degrees.
    scatter(zStarLocal(1,:)*(180/pi),zStarLocal(2,:)*(180/pi),'filled')
end
axis tight
h0=title('Stars Observed Every 5 Minutes for 20 Minutes');
h1=xlabel('Local Azimuth^\circ');
h2=ylabel('Local Elevation^\circ');
set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
set(h0,'FontSize',14,'FontWeight','bold','FontName','Times')
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
