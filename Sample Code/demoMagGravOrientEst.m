function demoMagGravOrientEst()
%%DEMOMAGGRAVOPRIENTEST  Demonstate how one can determine the relationshp
%                        between their local orientation with respect to
%                        the global WGS-84 coordinate system using a
%                        magnetic field measurement and a gravitational
%                        acceleration vector measurement. This can also be
%                        viewed as a general example of how one can use two
%                        sets of vectors to solve for a rotation matrix.
%
%Note that the accuracy of orientation estimation using magnetic and
%gravitational measurments is limited by the magnetic and gravitational
%models and can be quite bad at the geographic poles. Moreover, the
%function findTransParam, which is used to determine the relationship
%between different coordinate systems, becomes more subject to numerical
%precision errors as vectors become more co-linear.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

disp('Orientation Estimation Using Magnetic and Gravitational Measurements') 
disp('The observer is located in Hilo, Hawaii with zero ellipsoidal height.')

%The latitude and longitude of the Mauna Kea Observatory in Hilo, Hawaii.
phi=19.823*pi/180;%North latitude in radians.
lambda=-155.470*pi/180;%East longitude in radians.
height=0;
obsLoc=[phi;lambda;height];

disp('The true orientation of the observer is a rotation of 60 degrees')
disp('about the local ENU East (x)-axes 15 degrees about the rotated z axis.')
%We need to build up the true rotation matrix. First, we need to figure out
%the relationship between the local ENU axes and the global coordinate
%system.

%First, figure out how the ENU axes are related to the global WGS-84 axes.
%Get the local ENU axes. They are the local x-y-z axes.
uENU=getENUAxes(obsLoc);
%Get the rotation matrix to go from ECEF to ENU.
ECEF2ENURotMat=findTransParam(eye(3),uENU);

%Now, we need the rotation matrix from the local ENU axes to the global
%WGS-84 axes. It is assumed that the observer's true orientation is a
%rotation 60 degrees about the local ellipsoidal ENU x-axis-up followed by
%a rotation of 15 degrees about the rotated z-axis. So, we will now build
%the ENU2Local rotation matrix.
R1=axisAng2RotMat([1;0;0],60*pi/180);
R2=axisAng2RotMat([0;0;1],15*pi/180);
ENU2LocRotMat=R2*R1;

%Thus, the rotation matrix from global coordinates to local coordinates is
ECEF2LocRotMat=ENU2LocRotMat*ECEF2ENURotMat;

disp(' ')%Add a line break.
disp('1) Determining the gravitational acceleration and magnetic field at')
disp('the observer')
%The following vectors are determined in ECEF coordinates.

%Magnetic field coefficients for the EMM model.
[CMag,SMag,aMag,cMag]=getEMMCoeffs();
%Get the local magnetic field vector. It reqires the observer's location in
%spherical coordinates.
obsLocSpher=ellips2Sphere(obsLoc);
[~,negB]=spherHarmonicEval(CMag,SMag,obsLocSpher,aMag,cMag);
B=-negB;

%Get the zero-tide spherical harmonic coefficients.
[C,S,a,c]=getEGMGravCoeffs(2190,false);
[~,g]=spherHarmonicEval(C,S,obsLocSpher,a,c);%Acceleration due to gravity.

%The gravitational acceleration and magnetic field vectors as measured by
%the observer (ignoring noise) would be 
BLocal=ECEF2LocRotMat*B;
gLocal=ECEF2LocRotMat*g;

angDiffDeg=angBetweenVecs(B,g)*(180/pi);
disp(['The magnetic and gravitational acceleration vectors are ',num2str(angDiffDeg)])
disp('degrees apart.')

disp(' ')%Add a line break.
disp('2) Using the findTransparam function to determine orientation with')
disp('error-free measurements and NON-NORMALIZED vectors')
ECEF2LocRotMatEst=findTransParam([gLocal,BLocal],[g,B]);

xAxisEst=ECEF2LocRotMatEst*[0;0;1];
xAxisTrue=ECEF2LocRotMat*[0;0;1];
disp('The difference in the 3D rotation of the estimated and actual global')
disp(['x axes (due to numerical precision limitations) is ',num2str(angBetweenVecs(xAxisEst,xAxisTrue)*(180/pi)*60*60),' arcseconds'])

disp(' ')%Add a line break.
disp('3) Using the findTransparam function to determine orientation with')
disp('error-free measurements and NORMALIZED vectors')

gLocal=gLocal/norm(gLocal);
BLocal=BLocal/norm(BLocal);
g=g/norm(g);
B=B/norm(B);

ECEF2LocRotMatEst=findTransParam([gLocal,BLocal],[g,B]);

xAxisEst=ECEF2LocRotMatEst*[0;0;1];
xAxisTrue=ECEF2LocRotMat*[0;0;1];
disp('The difference in the 3D rotation of the estimated and actual global')
disp(['x axes (due to numerical precision limitations) is ',num2str(angBetweenVecs(xAxisEst,xAxisTrue)*(180/pi)*60*60),' arcseconds'])

disp('It can be seen that normalization of the vectors improves numeric')
disp('stability of deterministic orientation estimation')

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
