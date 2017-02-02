function [rICRS,rITRS]=solarBodyVec(Jul1,Jul2,timeCoordSys,Object,algorithm,obsState,stateCoordSys,deltaTTUT1,xpyp,dXdY)
%%SOLARBODYVEC   Get a vector from a near-Earth object or an object offset
%                from a solar system body or the solar system barycenter to
%                a given solar body, such as a planet or the sun. The
%                vector can be given in the International Celestial
%                Reference System (ICRS) or it can be given in the
%                international terrestrial reference system, an
%                Earth-Centered Earth-fixed coordinate system that
%                coincides with the WGS-84 standard. Refraction corrections
%                must be applied separately.
%
%INPUTS: Jul1,Jul2 Two parts of a Julian date given in the time system
%                  specified by the second argument. The units of the date
%                  are days. The full date is the sum of both terms. The
%                  date is broken into two parts to provide more bits of
%                  precision. It does not matter how the date is split.
%timeCoordSys  A parameter specifying the coordinate system of the Julian
%              date given by Jul1,Jul2. Possible values are
%              'UTC' Universal coordinated time.
%              'TT'  Terrestrial time.
%              'TDB' Barycentric dynamical time (can only be used if
%                    algorithm =0 or 1)
%              'TAI' International atomic time.
%              'TCG' Geocentric coordinate time.
%              'GPS' GPS time.
%               If the input is not given in TDB, then it must be converted
%               to TDB. In computing barycentric dynamical time, which is
%               necessary to use the ephemerides, the clock is assumed to
%               be collocated with the center of the Earth --thus, no
%               topocentric corrections are applied. This might result in
%               an often insignificant loss of precision.
%Object        The chosen celestial object to which one wants a vector. If
%              this parameter is omitted, then he sun is the chosen
%              celestial object. Object can take the values:
%              'SUN', 'MOON', 'SOLAR SYSTEM BARYCENTER', 'MERCURY',
%              'VENUS', 'MERCURY BARYCENTER', 'VENUS BARYCENTER',
%              'EARTH', 'EARTH-MOON BARYCENTER', 'MARS BARYCENTER',
%              'JUPITER BARYCENTER', 'SATURN BARYCENTER',
%              'URANUS BARYCENTER', 'NEPTUNE BARYCENTER',
%              'PLUTO BARYCENTER'
%              The barycenter terms refer to the center of mass of the
%              planet and its moons. With Mars, the barycenter of Mars is
%              very close to the center of mass of Mars, since the moons of
%              Mars are small. These are the values for which ephemerides
%              are provided in the DE430. If this parameter is omitted,
%              then 'SUN' is chosen by default.
%algorithm     An optional parameter that selects which algorithm to use.
%              The possibilities are:
%              0  (The default if this parameter is omitted) Use the DE430
%                 ephemerides and the functions from NASA's SPICE library
%                 to find the GCRS position of the object with corrections
%                 for light-time and aberration. Note that a special
%                 relativistic aberration correction is used rather than
%                 the Newtonian one that is part of the SPICE library.
%                 Gravitational bending and delay of light are not taken
%                 into account.
%              1  Use the DE430 ephemerides and the functions
%                 from NASA's SPICE library to find the GCRS position of
%                 the object without corrections for light-time and
%                 aberration.
%              2  Use the low-precision algorithms from the Astronomical
%                 Almanac to find the apparent GCRS position of the sun.
%                 An error will be returned if Object is anything other
%                 than 'SUN' or 'MOON'.
%obsState      The 3X1 or 6X1 location and velocity vector an observer in
%              Cartesian coordinates (in meters) in the coordinate system
%              selected by the next argument. The format is
%              [x;y;z;xdot;ydot;zdot]. If a 3X1 vector is passed, it is
%              assumed that the velocity is zero. If omitted or an empty
%              matrix is passed, a stationary observer at the origin of the
%              coordinate system specified by stateCoordSys is used.
%stateCoordSys The coordinate system that obsState is in. The possibilities
%              are
%              'GCRS' The geocentric celestial reference system.
%              'ITRS' The international terrestrial reference system.
%              If algorithm=1, then other values can be any of the values
%              listed above that are valid for Object as well as 'BCRS',
%              which is synonymous with 'SOLAR SYSTEM BARYCENTER'. When
%              using anything but 'ITRS', the coordinate system is aligned
%              with the International Celestial Reference System (ICRS),
%              but is centered at the given object. 'EARTH' is the same as
%              'GCRS'. If this parameter is omitted, then ITRS is chosen by
%              default.
%deltaTTUT1    An optional parameter specifying the difference between TT
%              and UT1 in seconds. Values of this parameter can be obtained
%              from
%http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%              or 
%http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%              If this parameter is omitted or an empty matrix is passed,
%              then the value provided by the function getEOP will be used
%              instead. This parameter as well as the following two are not
%              used in algorithm 2 if rITRS is not requested and
%              stateCoordSys='GCRS'.
%xpyp          xpyp=[xp;yp] are the polar motion coordinates in radians
%              including the effects of tides and librations. If this
%              parameter is omitted or if an empty matrix is passed, the
%              value from the function getEOP will be used.
%dXdY          dXdY=[dX;dY] are the celestial pole offsets with respect to
%              the IAU 2006/2000A precession/nutation model in radians If
%              this parameter is omitted or an empty matrix is passed, the
%              value from the function getEOP will be used.
%
%OUTPUTS: rICRS A vector in Cartesian ICRS coordinates whose origin is
%               aligned with the origin associated with stateCoordSys. If
%               stateCoordSys is 'GCRS', 'ITRS', or 'EARTH', then this is
%               in the GCRS coordinate system. The units are in meters from
%               the observer to the object with any light-time and
%               aberration corrections as requested. If algorithm 0 or 2 is
%               used, then this is a 3X1 position vector. If algorithm 1 is
%               used, then this is a 6X1 vector as it also contains the
%               velocity in meters per second.
%         rITRS A vector in Cartesian ITRS coordinates in meters from the
%               observer to the object with any light-time and aberration
%               corrections as requested.  If algorithm 0 or 2 is used,
%               then this is a 3X1 position vector. If algorithm 1 is used,
%               then this is a 6X1 vector as it also contains the velocity
%               in meters per second. If stateCoordSys is anything but
%               'GCRS', 'ITRS', or 'EARTH', then just an empty matrix is
%               returned for this output.
%
%If algorithm=0 or algorithm=1, the National Aeronautics and Space
%Administration's (NASA's) Navigation and Ancillary Information Facility's
%(NAIF's) SPICE toolkit's function cspice_spkezr is used with NASA's DE430
%ephemerides to determine the location of the object. Light-time and
%aberration corrections are applied. Those effects are defined in Chapter 7
%of [1]. Note that unlike the function iauAb in the IAU's Standards of
%Fundamental Astronomy library, no attempt is made to correct for the Sun's
%gravitational potential, which adds an error of up to 4 microarcseconds.
%Additionally, gravitational deflection from the sun and other planets is
%not taken into account.
%
%The low precision algorithm for the Sun and Moon that is used if
%algorithm=2 are taken from pages C5 and D22 of [2], where it is documented
%that the algorithm for the Sun is accurate to one arcminute for years from
%1950 to 2050 and the algorithm for the moon is  accurate to 0.3 degrees in
%ecliptic longitude and 0.2 degrees in ecliptic latitude for the years 1950
%to 2050.
%
%Note that although the function calls for the SPICE toolkit refer to the
%J2000.0 coordinate system, the DE430 ephemerides are NOT given in J2000.0
%coordinates. Rather, as documented in the comments to
%http://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf
%They are given aligned with the ICRS. Thus, no conversions between J2000.0
%and ICRS coordinates are necessary.
%
%REFERENCES:
%[1] S. E. Urban and K. P.Seidelmann, Eds.,Explanatory Supplement to the
%    Astronomical Almanac, 3rd ed. Mill Valley, CA: University Science
%    Books, 2013.
%[2] Department of Defense, Navy, Nautical Almanac Office, The Astronomical
%    Almanac for the Year 2014. Department of the Navy, 2013.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
   timeCoordSys='UTC'; 
end

if(nargin<4)
   Object='SUN';
end

if(nargin<5)
    algorithm=0;
end

if(nargin<6||isempty(obsState))
   obsState=[0;0;0;0;0;0]; 
end

if(size(obsState,1)==3)
    %Assume stationary of only a position is given.
    obsState=[obsState;zeros(3,1)];
end

if(nargin<7)
   stateCoordSys='ITRS'; 
end

if(nargin<8)
    deltaTTUT1=[];
end

if(nargin<9)
    xpyp=[];
end

if(nargin<10)
    dXdY=[];
end

if(strcmp(stateCoordSys,'BCRS'))
    stateCoordSys='SOLAR SYSTEM BARYCENTER';
end

if(algorithm~=0&&algorithm~=1&&~(strcmp(stateCoordSys,'ITRS')||strcmp(stateCoordSys,'GCRS')||strcmp(stateCoordSys,'EARTH')))
   error('State coordinate systems other than ITRS and GCRS can only be used with algorithms 0 and 1'); 
end

%If the functions from the MICE library should be used.
if(algorithm==0||algorithm==1)
    %Get the ephemerides.
    ScriptPath=mfilename('fullpath');
    ScriptFolder = fileparts(ScriptPath);
    
    %Load the ephemerides.
    cspice_furnsh([ScriptFolder, '/data/de430.bsp'])
    
    %The time needs to be in TDB. Perform a conversion, if it does not
    %depend on the chosen coordinate system. If given TT or UTC, convert to
    %TT and wait until the coordinate system for the location of the
    %observer is known.
    switch(timeCoordSys)
        case 'GPS'
            [Jul1,Jul2]=GPS2TT(Jul1,Jul2);
        case 'TCG'
            [Jul1,Jul2]=TCG2TT(Jul1,Jul2);
        case 'TT'
        case 'TDB'
        case 'TAI'
            [Jul1,Jul2]=TAI2TT(Jul1,Jul2);
        case 'UTC'
            [Jul1,Jul2]=UTC2TT(Jul1,Jul2);
        otherwise
            error('An unsupported coordinate system for the time is given.')
    end
    
    if(~strcmp(timeCoordSys,'TDB'))
        %Convert to TDB without topocentric corrections.
        [TDB1,TDB2]=TT2TDB(Jul1,Jul2,deltaTTUT1,[0;0;0]);
    else
        TDB1=Jul1;
        TDB2=Jul2;
    end

    %Convert to seconds past J2000.0 as a single double-precision floating
    %point number. The date of J2000.0 in TDB is 2451545.0. There are
    %precisely 86400 seconds in a Julian day.
    TDBSec=86400*((TDB1-2451545.0)+TDB2);
    
    switch(stateCoordSys)
        case 'GCRS'
            %The SPICE functions label the coordinate system as J2000.0,
            %but actually use the ICRS. Thus, no conversion is necessary.
            %obsPosGeo=ICRS2J2000F(obsState(1:3),Jul1,Jul2);
            %obsVelGeo=ICRS2J2000F(obsState(4:6),Jul1,Jul2);
            obsPosGeo=obsState(1:3);
            obsVelGeo=obsState(4:6);
            stateCoordSys='EARTH';
        case 'ITRS'
            %First, convert the position and velocity components of 
            %obsState into GCRS coordinates. The includes adding in the
            %velocity due to the Earth's rotation rate.
            obsGCRS=ITRS2GCRS(obsState,Jul1,Jul2,deltaTTUT1,xpyp,dXdY);

            %No rotation from GCRS into J2000.0 coordinates is necessary,
            %because despite how the SPICE toolkit is labeled, it is
            %actually using ICRS coordinates.
            %obsPosGeo=ICRS2J2000F(obsGCRS(1:3),Jul1,Jul2);
            %obsVelGeo=ICRS2J2000F(obsGCRS(4:6),Jul1,Jul2);
            obsPosGeo=obsGCRS(1:3);
            obsVelGeo=obsGCRS(4:6);
            stateCoordSys='EARTH';
        otherwise
            %Otherwise, assume that algorithm 1 is being used and the
            %origin for the state is something other than the Earth.
            obsPosGeo=obsState(1:3);
            obsVelGeo=obsState(4:6);
    end
    
    if(algorithm==1)
        %Do not apply light-time correction or aberration corrections.
        state = cspice_spkezr(Object, ...
                              TDBSec, ...
                              'J2000', ...
                              'NONE', ...
                              stateCoordSys);
        %Convert to meters and move to the observer.
        J2000Position=state(1:3)*1e3-obsPosGeo;
        J2000Velocity=state(4:6)*1e3;
    else
    %Given the state of the observer in ICRS coordinates (which the SPICE
    %toolkit mistakenly labels J2000.0 coordinates, we can compute the
    %apparent position taking into account light-time. Perform three
    %iterations of the light-time correction (ignoring gravitational
    %deflection of the light due to other planets or the Sun).
        c=Constants.speedOfLight;%Units of meters per second.
        LTEst=0;
    %Assume that it will converge in three iterations.
        for curIter=1:3
                state = cspice_spkezr(Object, ...
                                      TDBSec-LTEst, ...
                                      'J2000', ...
                                      'NONE', ...
                                      stateCoordSys);                    
    %The returned state is a vector of 3D position and velocity with units
    %of kilometers and kilometers per second. Convert to meters and meters
    %per second and then compute the time it would take light to travel
    %from the object to the observer
            J2000Position=state(1:3)*1e3-obsPosGeo;
            LTEst=norm(J2000Position)/c;
        end

     %Next, apply a special relativistic correction to the postion of the
     %object. For this, we need the instantaneous velocity vector of the
     %state coordinate system in the J2000.0 coordinate system when the
     %measurement is taken.
        SolarBodyState = cspice_spkezr('SOLAR SYSTEM BARYCENTER', ...
                                    TDBSec, ...
                                   'J2000', ...
                                   'NONE', ...
                                   stateCoordSys);
        %Convert to meters and meters per second.
        SolarBodyState=SolarBodyState*1e3;
        
        %The velocity of the observer in the J2000.0 dynamical coordinate
        %system is the velocity of the Earth plus the velocity of the
        %observer with respect to the Earth. The addition is performed in a
        %special relativistically correct manner.
        obsVelJ2000=relVecAdd(SolarBodyState(4:6),obsVelGeo);
        
        %Get the distance from the origin of the state coordinate system to
        %the Sun
        OrigSunState = cspice_spkezr('SUN',...
                                    TDBSec,...
                                   'J2000',...
                                   'NONE',...
                                   stateCoordSys);
        %The distance from the observer to the Sun.
        sunDist=norm(OrigSunState(1:3)*1e3-obsPosGeo);
        
        %Perform the aberration correction.
        J2000Position=aberrCorr(J2000Position,obsVelJ2000,sunDist);
        J2000Velocity=[];%Velocity is not returned.
    end
    %Had the SPICE toolkit actually provided coordinates in the J2000.0
    %dynamical system, it would be necessary to rotate them. However,
    %despite the labeling, the coordinates are in the ICRS.
    rICRS=[J2000Position;J2000Velocity];
    
    %Unload the ephemerides. Otherwise, repeated calls to this function
    %will load them again and again until the SPICE system exceeds the
    %limit of 1300 loaded kernels and has an error.
    cspice_unload([ScriptFolder, '/data/de430.bsp'])
else%Otherwise, use the low-precision algorithms from the Astronomical
    %Almanac for the sun and the moon.
    
    %The time should be converted to TT for these approximations.
    switch(timeCoordSys)
        case 'GPS'
            [Jul1,Jul2]=GPS2TT(Jul1,Jul2);
        case 'TCG'
            [Jul1,Jul2]=TCG2TT(Jul1,Jul2);
        case 'TAI'
            [Jul1,Jul2]=TAI2TT(Jul1,Jul2);
        case 'TT'
        case 'UTC'
            [Jul1,Jul2]=UTC2TT(Jul1,Jul2);
        otherwise
            error('An unsupported coordinate system for the time is given.')
    end
    
    switch(stateCoordSys)
        case 'GCRS'
            %The observer position and velocity vectors in GCRS just has to
            %be rotated to the orientation of the J2000.0 dynamical
            %coordinate system that is used in the SPICE functions.
            obsPosGCRS=obsState(1:3);
        case 'ITRS'
            %First, convert the position and velocity components of
            %obsState into GCRS coordinates. The includes adding in the
            %velocity due to the Earth's rotation rate.
            obsPosGCRS=ITRS2GCRS(obsState(1:3),Jul1,Jul2,deltaTTUT1,xpyp,dXdY);
        otherwise
            error('An unknown coordinate system is given for the state of the observer')
    end
    
    switch(Object)
        case 'SUN'
            n=(Jul1-2451545.0)+Jul2;
            %Mean longitude of the Sun in degrees, corrected for
            %aberration.
            L=mod(280.460+0.9856474*n,360);
            %Mean anomaly in degrees.
            g=mod(357.528+0.9856003*n,360);

            %The geocentric apparent ecliptic longitude in degrees adjusted
            %for aberration.
            lambda=L+1.915*sind(g)+0.020*sind(2*g);
            %The obliquity of the ecliptic in degrees.
            epsilon=23.439-0.0000004*n;

            %The distance of the sun from the Earth in astronomical units.
            R=1.00014-0.01671*cosd(g)-0.00014*cosd(2*g);

            %The value of rSun is in the true equinox and ecliptic of date
            %coordinate system denominated in astronomical units.
            rSun=[R*cosd(lambda);
                  R*cosd(epsilon)*sind(lambda);
                  R*sind(epsilon)*sind(lambda)];

            %Convert from astronomial units to meters.
            rSun=rSun*Constants.AstronomialUnit;

            %Transform rSun to the GCRS coordinate system.
            rICRS=TOD2GCRS(rSun,Jul1,Jul2);
            
            %Now, just subtract the position vectors to get the Sun at the
            %observer's location.
            rICRS=rICRS-obsPosGCRS;
        case 'MOON'
            %Julian centuries from J2000.0
            T=(Jul1-2451545)/36525;

            %The geocentric apparent ecliptic longitude in degrees.
            lambda=218.32+481267.881*T...
                   +6.29*sind(135.0+477198.87*T)-1.27*sind(259.3-413335.36*T)...
                   +0.66*sind(235.7+890534.22*T)+0.21*sind(269.9+954397.74*T)...
                   -0.19*sind(357.5+35999.05*T)-0.11*sind(186.5+966404.03*T);

            %The ecliptic latitude in degrees
            beta=5.13*sind(93.3+483202.02*T)+0.28*sind(228.2+960400.89*T)...
                -0.28*sind(318.3+6003.15*T)-0.17*sind(217.6-407332.21*T);

            piVal=0.9508+0.0518*cosd(135.0+477198.87*T)+0.0095*cosd(259.3-413335.36*T)...
                +0.0078*cosd(235.7+890534.22*T)+0.0028*cosd(269.9+954397.74*T);

            %The distance from the center of th Earth to the moon
            %denominated in Earth radii.
            r=1/sind(piVal);

            %Geodetic direction cosines
            l=cosd(beta)*cosd(lambda);
            m=0.9175*cosd(beta)*sind(lambda)-0.3978*sind(beta);
            n=0.3978*cosd(beta)*sind(lambda)+0.9175*sind(beta);

            %The value of rMoon is in the true equinox and ecliptic of date
            %coordinate system in units of Earth radii.
            rMoon=[r*l;
                   r*m;
                   r*n];

            %Convert from units of Earth radii to meters.
            rMoon=rMoon*Constants.EarthEqRadius;

            %Transform rMoon to the GCRS coordinate system.
            rICRS=TOD2GCRS(rMoon,Jul1,Jul2);
            
            %Now, just subtract the position vectors to get the Sun at the
            %observer's location.
            rICRS=rICRS-obsPosGCRS;
        otherwise
            error('Algorithm 2 can only be used with the sun and the moon')
    end
end

%If the location in ITRS coordinates is desired.
if(nargout>1)
    if(~(strcmp(stateCoordSys,'ITRS')||strcmp(stateCoordSys,'GCRS')||strcmp(stateCoordSys,'EARTH')))
        rITRS=[];
    else
        rITRS=GCRS2ITRS(rICRS,Jul1,Jul2,deltaTTUT1,xpyp,dXdY);
    end
end

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
