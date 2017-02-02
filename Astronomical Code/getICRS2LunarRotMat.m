function R=getICRS2LunarRotMat(Jul1,Jul2,timeCoordSys,moonCoordSys,deltaTTUT1)
%%GETICRS2LUNARROTMAT Get a rotation matrix to rotate a vector from the
%                International Celestial Refernce System (ICRS), whose axes
%                are aligned with those of the Barycentric Celestial
%                Reference System (BCRS) and the Geocentric Celestial
%                Reference System (GCRS) into the orientation of one
%                of two a lunar body-fixed coordinate systems, specified by
%                moonSys.
%
%INPUTS:Jul1,Jul2 Two parts of a Julian date given in the time system
%                  specified by the second argument. The units of the date
%                  are days. The full date is the sum of both terms. The
%                  date is broken into two parts to provide more bits of
%                  precision. It does not matter how the date is split.
%timeCoordSys  A parameter specifying the coordinate system of the Julian
%              date given by Jul1,Jul2. Possible values are
%              'UTC' Universal coordinated time.
%              'TT'  Terrestrial time.
%              'TDB' Barycentric dynamical time (can only be used if
%                    algorithm =0)
%              'TAI' International atomic time.
%              'TCG' Geocentric coordinate time.
%              'GPS' GPS time.
%               If the input is not given in TDB, then it must be converted
%               to TDB. In computing barycentric dynamical time, which is
%               necessary to use the ephemerides, the clock is assumed to
%               be collocated with the center of the Earth --thus, no
%               topocentric corrections are applied. This might result in a
%               minor loss of precision.
%moonCoordSys  A parameter indicating which of the two moon coordinate
%              systems to use. Possible values are
%              'MOON_PA' (the default if omitted or an empty matrix is
%                        passed) Use the principle axis coordinate system
%                        of the Moon. This is defined in terms of the
%                        Moon's inertia tensor. This is usually the
%                        coordinate system to use when dealing with
%                        spherical harmonic representations of the
%                        Moon's gravity. 
%              'MOON_ME' Use the Mean Earth/Polar Axis coordinate system
%                        for the Moon. This coordinate system is most
%                        appropriate when describing points on the surface
%                        of the Moon.
%%deltaTTUT1   An optional parameter that is only used if the time is not
%              given in TDB. Values of this parameter can be obtained
%              from
%http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%              or 
%http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%              If this parameter is omitted or an empty matrix is passed,
%              then the value provided by the function getEOP will be used
%              instead.
%               
%The two lunar coordinate systems are mentioned in [1].
%
%Routines from the the National Aeronautics and Space Administration's
%(NASA's) Navigation and Ancillary Information Facility's (NAIF's) SPICE
%toolkit are used to obtain the orientation of the Moon based on NASA's
%DE421 ephemerides. The DE42 ephemerides are used rather than the newer
%DE430 ephemerides, because the difference is not very big and no Moon
%orientation data with the DE430 ephemerides has come out for the SPICE
%toolkit yet.
%
%REFERENCES:
%[1] R. B. Roncoli, "Lunar constants and models document," Jet
%    Propulsion Laboratory, California Institute of Technology, Tech.
%    Rep. JPL D-32296, 23 Sep. 2005. [Online]. Available: 
%    http://www.hq.nasa.gov/alsj/lunar cmd 2005 jpl d32296.pdf
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3)
        timeCoordSys='TT';
    end

    if(nargin<4||isempty(moonCoordSys))
        moonCoordSys='MOON_PA';
    end
    
    %Load the data for the orientaiton of the moon. Two kernels must be
    %loaded.
    ScriptPath=mfilename('fullpath');
    ScriptFolder = fileparts(ScriptPath);
    cspice_furnsh([ScriptFolder, '/data/moon_pa_de421_1900-2050.bpc'])
    cspice_furnsh([ScriptFolder, '/data/moon_080317.tf'])
    
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

    if(~strcmp(timeCoordSys,'TDB')&&(nargin<5||isempty(deltaTTUT1)))
        [Jul1UTC,Jul2UTC]=TT2UTC(Jul1,Jul2);
        [~,~,~,deltaTTUT1]=getEOP(Jul1UTC,Jul2UTC);
    end
    
    %If the time has to be converted and no Earth orientation parameter is
    %provided, just use whatever the getEOP function returns.
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
    
    %Get the coordinate system
    R=cspice_pxform('J2000',moonCoordSys,TDBSec);
    
    %Unload the orientation data. Otherwise, repeated calls to this 
    %function will load them again and again until the SPICE system exceeds
    %the limit of 1300 loaded kernels and has an error.
    cspice_unload([ScriptFolder, '/data/moon_080317.tf'])
    cspice_unload([ScriptFolder, '/data/moon_pa_de421_1900-2050.bpc'])
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
