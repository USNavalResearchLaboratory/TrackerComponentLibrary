function DeltaR=solidTideShift(rStation,rSun,rMoon,Jul1,Jul2,addPermTide)
%%SOLIDTIDESHIFT Compute the vector offset of a point on the surface of the
%             Earth due to solid Earth tidal effects. These effects can be
%             over 31cm. If a Julian date is given, third-order solid Earth
%             tidal effects will be taken into account. Note that site
%             displacement due to ocean loading, which is not included in
%             the solid Earth tide offset, can also be quite significant,
%             with the IERS Conventions 2010 listing it up to 10cm.
%
%INPUTS: rStation A vector from the geocenter to the (non-offset) location
%                 of a station on the surface of the Earth in WGS-84 Earth-
%                 centered Earth-fixed (ECEF) coordinates in meters. This
%                 value is converted into spherical coordinates by this
%                 function, and only the angles in spherical coordinates,
%                 not the radius, are important for this algorithm.
%                 However, since locations on the ground are generally
%                 given in WGS-84 ellipsoidal coordinates, the ellipsoidal
%                 height of the ground plays a role in the conversion to
%                 spherical coordinates. However, variations in the
%                 ellipsoidal height provided on the order of 2km result in
%                 differences in the output DeltaR on the order of microns,
%                 so in practice, one could just set rStation to the
%                 Cartesian position corresponding to a particular
%                 ellipsoidal latitude and longitude with zero ellipsoidal
%                 height. Ideally, however, rStation is the Cartesian
%                 location of the land in a mean-tide model (if
%                 addPermTide=false), or in a tide-free model if
%                 addPermTide=true.
%            rSun A vector from the geocenter to the sun in ITRS
%                 coordinates in meters. This can be obtained using the
%                 readJPLEphem and GCRS2ITRS functions.
%           rMoon A vector from the geocenter to the moon in ITRS
%                 coordinates in meters. This can be obtained using the
%                 readJPLEphem and GCRS2ITRS functions.
%      Jul1, Jul2 Two parts of a Julian date given in terrestrial time
%                 (TT). The units of the date are days. The full date is
%                 the sum of both terms. The date is broken into two
%                 parts to provide more bits of precision. It does not
%                 matter how the date is split. These parameters must be
%                 given for the third-order tidal components to be taken
%                 into account.
%     addPermTide The deltaR values are for a conventional tide-free model.
%                 If addPermTide is true, then a constant (depending on
%                 latitude) permanent tide offset will be added to put the
%                 point into a mean-tide model. If this parameter is
%                 omitted, the default value is "false".
%
%OUTPUTS: DeltaR The offset of a point on the surface of the Earth due to
%                solid Earth tides. rStation+DeltaR is the location of the
%                station (point on the ground) at the terrestial time given
%                by Jul1 and Jul2 taking into account solid Earth tides.
%                Note that the addPermTide term must be consistent with the
%                coordinate system of rStation. If rStation already
%                includes the permanent tides (is in a mean-tide model),
%                then addPermTide should be false. Otherwise, addPermTide
%                should be true.
%
%The Sun and Moon cause the crust of the Earth to warp and points on the
%Earth in WGS-84 coordinates to move over time. The formulae for computing
%tidal shits of the crust are given in Section 7.1 of [1].
%
%Section 7.1.1 discusses solid Earth tides. While the second-order
%corrections are well documneted, the third order corrections can not be
%implemented directly from the standard. However, the standard specifies
%FORTRAN routines, whose algorithms do not resemble the work given in the
%standard. These Fortran routines, have been converted to Matlab and are
%called for the third-order corrections if Julian dates are given. The
%converted routines are separate functions at the bottom of this file.
%
%Note that third-order model has an implied set of ephemerides built into
%it that might not be perfectly consistent with whatever model was used to
%obtain rSun and rMoon.
%
%One time component is supposed to be in Julian centuries since J2000.0. As
%J2000.0 is defined in TDB and the time is given in TT, the TDB2TT function
%is used. The extra parameters needed for the TDB2TT function are not
%requested as inputs to this function, since it is assumed that the
%difference between the tides is too small to matter for the model.
%
%To test this function, one can use the same values that are given in the
%IERS's implementation, DEHANTTIDEINEL.F Those are
% rStation = [4075578.385;931852.890;4801570.154];
% rSun=[137859926952.015;54228127881.4350;23509422341.6960];
% rMoon=[-179996231.920342;-312468450.131567;-169288918.592160];
% addPermTide=false;
% [Jul1,Jul2]=Cal2TT(2009,4,13,0,0,0);
% DeltaR=solidTideShift(rStation,rSun,rMoon,addPermTide,Jul1,Jul2);
%
%The value of DeltaR truncated to 9 places should be about
%DeltaR = 0.0770042035
%         0.0630405632
%         0.0551656815
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The convention for indexing arrays is that 1 refers to a value with
%respect to lunar parameters and 2 refers to a value with respect to solar
%parameters.

if(nargin<6||isempty(addPermTide))
    addPermTide=false;
end

%%%%%%General constants that are used%%%%%%
Re=Constants.EarthEqRadius;%meters
GMGMe(1)=Constants.MoonEarthMassRatio;
GMGMe(2)=Constants.SunEarthMassRatio;

%%%%%%Convert the input parameters into different coordinate systems%%%%%%

%Obtain the spherical coordinates of the station.
pointSpherStation=Cart2Sphere(rStation);
lambda=pointSpherStation(2);%Geocentric longitude
phi=pointSpherStation(3);%Geocentric latitude.
r=pointSpherStation(1);
rHat=rStation/r;%A unit vector in the direction of the station.

%Get East-North-Up Axes using a local SPHERICAL Earth model. A flattening
%factor of zero makes the Earth flat.
uENU=getENUAxes([phi;lambda;0],false,r,0);
uE=uENU(:,1);%Spherical East unit vector
uN=uENU(:,2);%Spherical North unit vector.
uU=uENU(:,3);%Spherical radial *up) unit vector

%Extract the x,y, and z coordinates of the moon.
Xj(1)=rMoon(1);
Yj(1)=rMoon(2);
Zj(1)=rMoon(3);
%Obtain the spherical coordinates of the moon.
pointSpherMoon=Cart2Sphere(rMoon);
Lambdaj(1)=pointSpherMoon(2);
Phij(1)=pointSpherMoon(3);
Rj(1)=pointSpherMoon(1);
RjHat(:,1)=rMoon/Rj(1);%A unit vector in the direction of the moon.

%Extract the x,y, and z coordinates of the sun.
Xj(2)=rSun(1);
Yj(2)=rSun(2);
Zj(2)=rSun(3);
%Obtain the spherical coordinates of the sun.
pointSpherSun=Cart2Sphere(rSun);
Lambdaj(2)=pointSpherSun(2);
Phij(2)=pointSpherSun(3);
Rj(2)=pointSpherSun(1);
RjHat(:,2)=rSun/Rj(2);%A unit vector in the direction of the sun.

%%%%%%Determine the solid Earth tidal displacement.%%%%%%
%This is a three step process. The first two steps are described in the
%table on page 103 of the 2010 IERS conventions.

%%%Step 1: time-domain contributions
%%First, the in-phase contribution.

%Determine the second degree nomial Love (h) and Shida (l) numbers using
%the latitude correction in the second paragraph on page 105.
h0=0.6078;
h2=-0.0006;
l0=0.0847;
l2=0.0002;

h2=h0+h2*(3*sin(phi)^2-1)/2;
l2=l0+l2*(3*sin(phi)^2-1)/2;

%Equation 7.5 for the degree 2 tidal offset
DeltaR=0;
for n=1:2
    RjHatDotrHat=dot(RjHat(:,n),rHat);
    DeltaR=DeltaR+GMGMe(n)*(Re^4/Rj(n)^3)*(h2*rHat*(3*RjHatDotrHat^2-1)/2+3*l2*RjHatDotrHat*(RjHat(:,n)-RjHatDotrHat*rHat));
end

%Equation 7.6 for the degree 3 tidal offset
h3=0.292;
l3=0.015;
for n=1:2
    RjHatDotrHat=dot(RjHat(:,n),rHat);
    DeltaR=DeltaR+GMGMe(n)*(Re^5/Rj(n)^4)*(h3*rHat*((5/2)*RjHatDotrHat^3-(3/2)*RjHatDotrHat)+l3*((15/2)*RjHatDotrHat^2-3/2)*(RjHat(:,n)-RjHatDotrHat*rHat));
end

%Equation 7.8 for the diurnal latitude dependence of the Love numbers.
l1Diurn=0.0012;
deltaT=0;
for n=1:2
    P12Cos=3*Xj(n)*Zj(n)/Rj(n)^2;%Equation 7.7b
    P12Sin=3*Yj(n)*Zj(n)/Rj(n)^2;%Equation 7.7b
    
    %This is a term in the sum of Equation 7.8. However, it does not look
    %like Equation 7.8, since angle difference formulae had to be used to
    %get terms suitable for the values in Equation 7.7b to be used.
    %Specifically,
    %sin(lambda-Lambdaj(n))=sin(lambda)*cos(Lambdaj(n))-cos(lambda)*sin(Lambdaj(n));
    %cos(lambda-Lambdaj(n))=cos(lambda)*cos(Lambdaj(n))+sin(lambda)*sin(Lambdaj(n));
    deltaT=deltaT+GMGMe(n)*(Re^4/Rj(n)^3)*(uN*sin(phi)*(cos(lambda)*P12Cos+sin(lambda)*P12Sin)-uE*cos(2*phi)*(sin(lambda)*P12Cos-cos(lambda)*P12Sin));
end
deltaT=-l1Diurn*sin(phi)*deltaT;

DeltaR=DeltaR+deltaT;

%Equation 7.9 for the semi-diurnal latitude dependence of the Love numbers.
l1SemiDiurn=0.0024;
deltaT=0;
for n=1:2
    P22Cos=(3/Rj(n)^2)*(Xj(n)^2-Yj(n)^2);%Equation 7.7c
    P22Sin=(6/Rj(n)^2)*Xj(n)*Yj(n);%Equation 7.7c
    
    %This is a term in the sum of Equation 7.9. However, it does not look
    %like Equation 7.9, since angle difference formulae had to be used to
    %get terms suitable for the values in Equation 7.7c to be used.
    %Specifically,
    %sin(2*(lambda-Lambdaj(n)))=sin(2*lambda)*cos(2*Lambdaj(n))-cos(2*lambda)*sin(2*Lambdaj(n));
    %cos(2*(lambda-Lambdaj(n)))=cos(2*lambda)*cos(2*Lambdaj(n))+sin(2*lambda)*sin(2*Lambdaj(n));
    deltaT=deltaT+GMGMe(n)*(Re^4/Rj(n)^3)*(uN*(cos(2*lambda)*P22Cos+sin(2*lambda)*P22Sin)+uE*sin(phi)*(sin(2*lambda)*P22Cos-cos(2*lambda)*P22Sin));
end
deltaT=-(1/2)*l1SemiDiurn*sin(phi)*cos(phi)*deltaT;

DeltaR=DeltaR+deltaT;

%%Next, the out-of-phase contributions for the second degree terms

%First, the diurnal tides in Equations 7.10a and 7.10b
hIDiurn=-0.0025;
lIDirun=-0.0007;

deltaR=0;
deltaT=0;
for n=1:2
    deltaR=deltaR+GMGMe(n)*(Re^4/Rj(n)^3)*sin(2*Phij(n))*sin(2*phi)*sin(lambda-Lambdaj(n));
    deltaT=deltaT+GMGMe(n)*(Re^4/Rj(n)^3)*sin(2*Phij(n))*(cos(2*phi)*sin(lambda-Lambdaj(n))*uN+sin(phi)*cos(lambda-Lambdaj(n))*uE);
end
deltaR=-(3/4)*hIDiurn*deltaR;
deltaT=-(3/2)*lIDirun*deltaT;

DeltaR=DeltaR+deltaR*uU+deltaT;

%Next, the semidiurnal tides in Equations 7.11a and 7.11b.
hISemiDiurn=-0.0022;
lISemiDiurn=-0.0007;

deltaR=0;
deltaT=0;
for n=1:2
    deltaR=deltaR+GMGMe(n)*(Re^4/Rj(n)^3)*cos(Phij(n))^2*cos(phi)^2*sin(2*(lambda-Lambdaj(n)));
    deltaT=deltaT+GMGMe(n)*(Re^4/Rj(n)^3)*cos(Phij(n))^2*(sin(2*phi)*sin(2*(lambda-Lambdaj(n)))*uN-2*cos(phi)*cos(2*(lambda-Lambdaj(n)))*uE);
end
deltaR=-(3/4)*hISemiDiurn*deltaR;
deltaT=(3/4)*lISemiDiurn*deltaT;

DeltaR=DeltaR+deltaR*uU+deltaT;

%%%Step 2: frequency-domain contributions

%Only apply the corrections if a Julian date has been provided.
if(nargin>4)
    %%First, the contributions for the diurnal band      
      %Convert Terrestrial Time to Julian centuries since J2000.0. A Julian
      %century is defined to have 36525 days in it. The offset of 2451545.0
      %days is the Julian date at J2000.0 in TDB. In TT, a difference of a
      %few milliseconds exists.
      [TDB1,TDB2]=TT2TDB(Jul1,Jul2);
      T=((TDB1-2451545.0)+TDB2)/36525;

      %Convert the Julian date to UTC and determine the fractional number
      %of hours passed in the day assuming precisely 24 hours in a say.
      [Jul1,Jul2]=TT2UTC(Jul1,Jul2);
      [~,~,~,dayFrac]=UTC2Cal(Jul1,Jul2,true);
      FHR=24*dayFrac;

    DeltaR=DeltaR+diurBandCorr(rStation,FHR,T);

    %%Next, the contributions from the long-period band
    DeltaR=DeltaR+longBandCorr(rStation,T);
end

if(addPermTide~=false)
    P2=(3*sin(phi)^2-1)/2;
    %Equation 7.14a in the IERS conventions
    deltaR=(-0.1206+0.0001*P2)*P2;
    %Equation 7.14b in the IERS conventions
    deltaT=(-0.0252-0.0001*P2)*sin(2*phi);
    DeltaR=DeltaR+deltaR*uU+deltaT*uN;
end
end


function XCORSTA=diurBandCorr(XSTA,FHR,T)
%%DIURBANDCORR  This subroutine is a Matlab translation of the subroutine
%               STEP2DIU that is provided by the IERS at
%               ftp://tai.bipm.org/iers/convupdt/chapter7/
%               The algorithm is not fully documented by the IERS 2010
%               conventions. For example, most of the constants are not
%               listed anywhere in the documentation.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%
%SUBROUTINE STEP2DIU (XSTA,FHR,T,XCORSTA)  
%
%  - - - - - - - - - - -
%   S T E P 2 D I U
%  - - - - - - - - - - -
%
%  This routine is part of the International Earth Rotation and
%  Reference Systems Service (IERS) Conventions software collection.
%
%  This subroutine gives the in-phase and out-of-phase corrections
%  induced by mantle anelasticity in the diurnal band. 
%
%  In general, Class 1, 2, and 3 models represent physical effects that
%  act on geodetic parameters while canonical models provide lower-level
%  representations or basic computations that are used by Class 1, 2, or
%  3 models.
% 
%  Status: Class 1
%
%     Class 1 models are those recommended to be used a priori in the
%     reduction of raw space geodetic data in order to determine
%     geodetic parameter estimates.
%     Class 2 models are those that eliminate an observational
%     singularity and are purely conventional in nature.
%     Class 3 models are those that are not required as either Class
%     1 or 2.
%     Canonical models are accepted as is and cannot be classified as a
%     Class 1, 2, or 3 model.
%
%  Given:
%     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
%     FHR           d      Fractional hours in the day (Note 2)
%     T             d      Centuries since J2000
%
%  Returned:
%     XCORSTA       d(3)   In phase and out of phase station corrections
%                          for diurnal band (Note 4)
%
%  Notes:
%
%  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
%     expressed in meters. 
%  
%  2) The fractional hours in the day is computed as the hour + minutes/60.0
%     + sec/3600.0.  The unit is expressed in Universal Time (UT).
%
%  4) All coordinates are expressed in meters.
%
%  Test case:
%     given input: XSTA(1) = 4075578.385D0 meters
%                  XSTA(2) =  931852.890D0 meters
%                  XSTA(3) = 4801570.154D0 meters 
%                  FHR     = 0.00D0 hours
%                  T       = 0.1059411362080767D0 Julian centuries
%                  
%     expected output:  XCORSTA(1) = 0.4193085327321284701D-02 meters
%                       XCORSTA(2) = 0.1456681241014607395D-02 meters
%                       XCORSTA(3) = 0.5123366597450316508D-02 meters
%
%  References:
%
%     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
%     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
%
%     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
%     IERS Technical Note No. 36, BKG (2010)
%
%  Revisions:
%  1996 March    23 V. Dehant      Original code
%  2009 July     31 B.E. Stetzler  Initial standardization of code 
%  2009 August   06 B.E. Stetzler  Provided a test case
%  2009 August   06 B.E. Stetzler  Capitalized all variables for 
%                                  Fortran 77 compatibility
%  2010 October  20 B.E. Stetzler  Input T corrected to be number of
%                                  centuries since J2000
%-----------------------------------------------------------------------
      
      D2PI = 6.283185307179586476925287;
      
      DATDI=[-3,    0,  2,   0,  0, -0.01,   0,      0,      0;
             -3,    2,  0,   0,  0, -0.01,   0,      0,      0;
             -2,    0,  1,  -1,  0, -0.02,   0,      0,      0;
             -2,    0,  1,   0,  0, -0.08,   0,     -0.01,   0.01;
             -2,    2, -1,   0,  0, -0.02,   0,      0,      0;
             -1,    0,  0,  -1,  0, -0.10,   0,      0,      0;
             -1,    0,  0,   0,  0, -0.51,   0,     -0.02,   0.03;
             -1,    2,  0,   0,  0,  0.01,   0,      0,      0;
              0,   -2,  1,   0,  0,  0.01,   0,      0,      0;
              0,    0, -1,   0,  0,  0.02,   0,      0,      0;
              0,    0,  1,   0,  0,  0.06,   0,      0,      0;
              0,    0,  1,   1,  0,  0.01,   0,      0,      0;
              0,    2, -1,   0,  0,  0.01,   0,      0,      0;
              1,   -3,  0,   0,  1, -0.06,   0,      0,      0;
              1,   -2,  0,  -1,  0,  0.01,   0,      0,      0;
              1,   -2,  0,   0,  0, -1.23,  -0.07,   0.06,   0.01;
              1,   -1,  0,   0, -1,  0.02,   0,      0,      0;
              1,   -1,  0,   0,  1,  0.04,   0,      0,      0;
              1,    0,  0,  -1,  0, -0.22,   0.01,   0.01,   0;
              1,    0,  0,   0,  0, 12.00,  -0.80,  -0.67,  -0.03;
              1,    0,  0,   1,  0,  1.73,  -0.12,  -0.10,   0;
              1,    0,  0,   2,  0, -0.04,   0,      0,      0;
              1,    1,  0,   0, -1, -0.50,  -0.01,   0.03,   0;
              1,    1,  0,   0,  1,  0.01,   0,      0,      0;
              0,    1,  0,   1, -1, -0.01,   0,      0,      0;
              1,    2, -2,   0,  0, -0.01,   0,      0,      0;
              1,    2,  0,   0,  0, -0.11,   0.01,   0.01,   0;
              2,   -2,  1,   0,  0, -0.01,   0,      0,      0;
              2,    0, -1,   0,  0, -0.02,   0,      0,      0;
              3,    0,  0,   0,  0,  0,      0,      0,      0;
              3,    0,  0,   1,  0,  0,      0,      0,      0];
      
    DEG2RAD = D2PI/360;

%  Compute the phase angles in degrees.
    S = 218.31664563+(481267.88194+(-0.0014663889+(0.00000185139)*T)*T)*T;

    TAU = FHR*15+280.4606184+(36000.7700536+(0.00038793+(-0.0000000258)*T)*T)*T+(-S);

    PR = (1.396971278+(0.000308889+(0.000000021+(0.000000007)*T)*T)*T)*T;

    S = S + PR;

    H = 280.46645+(36000.7697489+(0.00030322222+(0.000000020+(-0.00000000654)*T)*T)*T)*T;

    P = 83.35324312+(4069.01363525+(-0.01032172222+(-0.0000124991+(0.00000005263)*T)*T)*T)*T;

    ZNS = 234.95544499+(1934.13626197+(-0.00207561111+(-0.00000213944+(0.00000001650)*T)*T)*T)*T;

    PS = 282.93734098+(1.71945766667+(0.00045688889+(-0.00000001778+(-0.00000000334)*T)*T)*T)*T;

% Reduce angles to between the range 0 and 360.
    S =  mod(S,360);
    TAU = mod(TAU,360);
    H =  mod(H,360);
    P =  mod(P,360);
    ZNS = mod(ZNS,360);
    PS = mod(PS,360);

    RSTA = norm(XSTA); 
    SINPHI = XSTA(3)/RSTA;
    COSPHI = norm(XSTA(1:2))/RSTA;

    COSLA = XSTA(1)/COSPHI/RSTA;
    SINLA = XSTA(2)/COSPHI/RSTA;
    ZLA = atan2(XSTA(2),XSTA(1));
 
% Initialize.
    XCORSTA=zeros(3,1);

    for J=1:31
        % Convert from degrees to radians.
        THETAF=(TAU+DATDI(J,1)*S+DATDI(J,2)*H+DATDI(J,3)*P+DATDI(J,4)*ZNS+DATDI(J,5)*PS)*DEG2RAD;

        DR=DATDI(J,6)*2*SINPHI*COSPHI*sin(THETAF+ZLA)+DATDI(J,7)*2*SINPHI*COSPHI*cos(THETAF+ZLA);

        DN=DATDI(J,8)*(COSPHI^2-SINPHI^2)*sin(THETAF+ZLA)+DATDI(J,9)*(COSPHI^2-SINPHI^2)*cos(THETAF+ZLA);
        %      DE=DATDI(8,J)*SINPHI*COS(THETAF+ZLA)+
        %     Modified 20 June 2007

        DE=DATDI(J,8)*SINPHI*cos(THETAF+ZLA)-DATDI(J,9)*SINPHI*sin(THETAF+ZLA);

        XCORSTA(1)=XCORSTA(1)+DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA;
        XCORSTA(2)=XCORSTA(2)+DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA;  
        XCORSTA(3)=XCORSTA(3)+DR*SINPHI+DN*COSPHI;
    end

    XCORSTA=XCORSTA/1000;

%  Finished.

%+----------------------------------------------------------------------
%
%  Copyright (C) 2008
%  IERS Conventions Center
%
%  ==================================
%  IERS Conventions Software License
%  ==================================
%
%  NOTICE TO USER:
%
%  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
%  WHICH APPLY TO ITS USE.
%
%  1. The Software is provided by the IERS Conventions Center ("the
%     Center").
%
%  2. Permission is granted to anyone to use the Software for any
%     purpose, including commercial applications, free of charge,
%     subject to the conditions and restrictions listed below.
%
%  3. You (the user) may adapt the Software and its algorithms for your
%     own purposes and you may distribute the resulting "derived work"
%     to others, provided that the derived work complies with the
%     following requirements:
%
%     a) Your work shall be clearly identified so that it cannot be
%        mistaken for IERS Conventions software and that it has been
%        neither distributed by nor endorsed by the Center.
%
%     b) Your work (including source code) must contain descriptions of
%        how the derived work is based upon and/or differs from the
%        original Software.
%
%     c) The name(s) of all modified routine(s) that you distribute
%        shall be changed.
% 
%     d) The origin of the IERS Conventions components of your derived
%        work must not be misrepresented; you must not claim that you
%        wrote the original Software.
%
%     e) The source code must be included for all routine(s) that you
%        distribute.  This notice must be reproduced intact in any
%        source distribution. 
%
%  4. In any published work produced by the user and which includes
%     results achieved by using the Software, you shall acknowledge
%     that the Software was used in obtaining those results.
%
%  5. The Software is provided to the user "as is" and the Center makes
%     no warranty as to its use or performance.   The Center does not
%     and cannot warrant the performance or results which the user may
%     obtain by using the Software.  The Center makes no warranties,
%     express or implied, as to non-infringement of third party rights,
%     merchantability, or fitness for any particular purpose.  In no
%     event will the Center be liable to the user for any consequential,
%     incidental, or special damages, including any lost profits or lost
%     savings, even if a Center representative has been advised of such
%     damages, or for any claim by any third party.
%
%  Correspondence concerning IERS Conventions software should be
%  addressed as follows:
%
%                     Gerard Petit
%     Internet email: gpetit[at]bipm.org
%     Postal address: IERS Conventions Center
%                     Time, frequency and gravimetry section, BIPM
%                     Pavillon de Breteuil
%                     92312 Sevres  FRANCE
%
%     or
%
%                     Brian Luzum
%     Internet email: brian.luzum[at]usno.navy.mil
%     Postal address: IERS Conventions Center
%                     Earth Orientation Department
%                     3450 Massachusetts Ave, NW
%                     Washington, DC 20392
%
%
%-----------------------------------------------------------------------
end  


function XCORSTA=longBandCorr(XSTA,T)
%%DIURBANDCORR  This subroutine is a Matlab translation of the subroutine
%               STEP2LON that is provided by the IERS at
%               ftp://tai.bipm.org/iers/convupdt/chapter7/
%               The algorithm is not fully documented by the IERS 2010
%               conventions. For example, most of the constants are not
%               listed anywhere in the documentation.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%
%SUBROUTINE STEP2LON (XSTA,T,XCORSTA)
%
%  - - - - - - - - - - -
%   S T E P 2 L O N
%  - - - - - - - - - - -
%
%  This routine is part of the International Earth Rotation and
%  Reference Systems Service (IERS) Conventions software collection.
%
%  This subroutine gives the in-phase and out-of-phase corrections
%  induced by mantle anelasticity in the long period band. 
%
%  In general, Class 1, 2, and 3 models represent physical effects that
%  act on geodetic parameters while canonical models provide lower-level
%  representations or basic computations that are used by Class 1, 2, or
%  3 models.
% 
%  Status: Class 1
%
%     Class 1 models are those recommended to be used a priori in the
%     reduction of raw space geodetic data in order to determine
%     geodetic parameter estimates.
%     Class 2 models are those that eliminate an observational
%     singularity and are purely conventional in nature.
%     Class 3 models are those that are not required as either Class
%     1 or 2.
%     Canonical models are accepted as is and cannot be classified as a
%     Class 1, 2, or 3 model.
%
%  Given:
%     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
%     T             d      Centuries since J2000
%
%  Returned:
%     XCORSTA       d(3)   In phase and out of phase station corrections
%                          for diurnal band (Note 2)
%
%  Notes:
%
%  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
%     expressed in meters. 
%  
%  2) All coordinates are expressed in meters.
%
%  Test case:
%     given input: XSTA(1) = 4075578.385D0 meters
%                  XSTA(2) =  931852.890D0 meters
%                  XSTA(3) = 4801570.154D0 meters 
%                  T       = 0.1059411362080767D0 Julian centuries
%                  
%     expected output:  XCORSTA(1) = -0.9780962849562107762D-04 meters
%                       XCORSTA(2) = -0.2236349699932734273D-04 meters
%                       XCORSTA(3) =  0.3561945821351565926D-03 meters
%
%  References:
%
%     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
%     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
%
%     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
%     IERS Technical Note No. 36, BKG (2010)
%
%  Revisions:
%  1996 March    23 V. Dehant      Original code
%  2009 August   07 B.E. Stetzler  Initial standardization of code
%                                  and found unnecessary variables tau
%                                  and fhr 
%  2009 August   07 B.E. Stetzler  Provided a test case
%  2009 August   07 B.E. Stetzler  Capitalized all variables for 
%                                  Fortran 77 compatibility
%  2010 October  20 B.E. Stetzler  Input T corrected to be number of 
%                                  centuries since J2000
%-----------------------------------------------------------------------

    D2PI = 6.283185307179586476925287;

    DATDI=[0, 0, 0, 1, 0,   0.47, 0.23, 0.16, 0.07;
         0, 2, 0, 0, 0,  -0.20,-0.12,-0.11,-0.05;
         1, 0,-1, 0, 0,  -0.11,-0.08,-0.09,-0.04;
         2, 0, 0, 0, 0,  -0.13,-0.11,-0.15,-0.07;
         2, 0, 0, 1, 0,  -0.05,-0.05,-0.06,-0.03];


    DEG2RAD = D2PI/360;

%  Compute the phase angles in degrees.
    S = 218.31664563+(481267.88194+(-0.0014663889+(0.00000185139)*T)*T)*T;

    PR = (1.396971278+(0.000308889+(0.000000021+(0.000000007)*T)*T)*T)*T;

    S = S + PR;

    H = 280.46645+(36000.7697489+(0.00030322222+(0.000000020+(-0.00000000654)*T)*T)*T)*T; 

    P = 83.35324312+(4069.01363525+(-0.01032172222+(-0.0000124991+(0.00000005263)*T)*T)*T)*T;

    ZNS = 234.95544499+(1934.13626197+(-0.00207561111+(-0.00000213944+(0.00000001650)*T)*T)*T)*T;

    PS = 282.93734098+(1.71945766667+(0.00045688889+(-0.00000001778+(-0.00000000334)*T)*T)*T)*T;

    RSTA=norm(XSTA);
    SINPHI=XSTA(3)/RSTA;
    COSPHI=norm(XSTA(1:2))/RSTA;
    
    COSLA=XSTA(1)/COSPHI/RSTA;
    SINLA=XSTA(2)/COSPHI/RSTA;

% Reduce angles to between the range 0 and 360.
    S =  mod(S,360);
    %      TAU = DMOD(TAU,360D0)
    H =  mod(H,360);
    P =  mod(P,360);
    ZNS = mod(ZNS,360);
    PS = mod(PS,360);

    DR_TOT = 0;
    DN_TOT = 0;

    XCORSTA=zeros(3,1);
    
    for J=1:5
        THETAF=(DATDI(J,1)*S+DATDI(J,2)*H+DATDI(J,3)*P+DATDI(J,4)*ZNS+DATDI(J,5)*PS)*DEG2RAD;

        DR=DATDI(J,6)*(3D0*SINPHI^2-1)/2*cos(THETAF)+DATDI(J,8)*(3D0*SINPHI^2-1)/2*sin(THETAF);

        DN=DATDI(J,7)*(COSPHI*SINPHI*2)*cos(THETAF)+DATDI(J,9)*(COSPHI*SINPHI*2)*sin(THETAF);

        DE = 0;
        DR_TOT = DR_TOT+DR;
        DN_TOT = DN_TOT+DN;

        XCORSTA(1)=XCORSTA(1)+DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA; 
        XCORSTA(2)=XCORSTA(2)+DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA;
        XCORSTA(3)=XCORSTA(3)+DR*SINPHI+DN*COSPHI;
    end   

    XCORSTA=XCORSTA/1000;

%  Finished.

%+----------------------------------------------------------------------
%
%  Copyright (C) 2008
%  IERS Conventions Center
%
%  ==================================
%  IERS Conventions Software License
%  ==================================
%
%  NOTICE TO USER:
%
%  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
%  WHICH APPLY TO ITS USE.
%
%  1. The Software is provided by the IERS Conventions Center ("the
%     Center").
%
%  2. Permission is granted to anyone to use the Software for any
%     purpose, including commercial applications, free of charge,
%     subject to the conditions and restrictions listed below.
%
%  3. You (the user) may adapt the Software and its algorithms for your
%     own purposes and you may distribute the resulting "derived work"
%     to others, provided that the derived work complies with the
%     following requirements:
%
%     a) Your work shall be clearly identified so that it cannot be
%        mistaken for IERS Conventions software and that it has been
%        neither distributed by nor endorsed by the Center.
%
%     b) Your work (including source code) must contain descriptions of
%        how the derived work is based upon and/or differs from the
%        original Software.
%
%     c) The name(s) of all modified routine(s) that you distribute
%        shall be changed.
% 
%     d) The origin of the IERS Conventions components of your derived
%        work must not be misrepresented; you must not claim that you
%        wrote the original Software.
%
%     e) The source code must be included for all routine(s) that you
%        distribute.  This notice must be reproduced intact in any
%        source distribution. 
%
%  4. In any published work produced by the user and which includes
%     results achieved by using the Software, you shall acknowledge
%     that the Software was used in obtaining those results.
%
%  5. The Software is provided to the user "as is" and the Center makes
%     no warranty as to its use or performance.   The Center does not
%     and cannot warrant the performance or results which the user may
%     obtain by using the Software.  The Center makes no warranties,
%     express or implied, as to non-infringement of third party rights,
%     merchantability, or fitness for any particular purpose.  In no
%     event will the Center be liable to the user for any consequential,
%     incidental, or special damages, including any lost profits or lost
%     savings, even if a Center representative has been advised of such
%     damages, or for any claim by any third party.
%
%  Correspondence concerning IERS Conventions software should be
%  addressed as follows:
%
%                     Gerard Petit
%     Internet email: gpetit[at]bipm.org
%     Postal address: IERS Conventions Center
%                     Time, frequency and gravimetry section, BIPM
%                     Pavillon de Breteuil
%                     92312 Sevres  FRANCE
%
%     or
%
%                     Brian Luzum
%     Internet email: brian.luzum[at]usno.navy.mil
%     Postal address: IERS Conventions Center
%                     Earth Orientation Department
%                     3450 Massachusetts Ave, NW
%                     Washington, DC 20392
%
%
%-----------------------------------------------------------------------
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

