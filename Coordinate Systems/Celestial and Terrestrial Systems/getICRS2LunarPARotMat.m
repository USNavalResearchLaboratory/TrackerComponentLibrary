function R=getICRS2LunarPARotMat(Jul1,Jul2,timeCoordSys,deltaTTUT1)
%%GETICRS2LUNARPAROTMAT Get a rotation matrix to rotate a vector from the
%                International Celestial Refernce System (ICRS), whose axes
%                are aligned with those of the Barycentric Celestial
%                Reference System (BCRS) and the Geocentric Celestial
%                Reference System (GCRS) into the principle axis coordinate
%                system of the Moon. This is usually the coordinate system
%                to use when dealing with spherical harmonic
%                representations of the Moon's gravity. 
%
%INPUTS:Jul1,Jul2 Two parts of a Julian date given in the time system
%               specified by the second argument. The units of the date are
%               days. The full date is the sum of both terms. The date is
%               broken into two parts to provide more bits of precision. It
%               does not matter how the date is split.
% timeCoordSys A parameter specifying the coordinate system of the Julian
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
%   deltaTTUT1 An optional parameter that is only used if the time is not
%              given in TDB. Values of this parameter can be obtained
%              from
%              http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%              or 
%              http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%              If this parameter is omitted or an empty matrix is passed,
%              then the value provided by the function getEOP will be used
%              instead.
%               
%Two lunar coordinate systems are mentioned in [1].
%
%REFERENCES:
%[1] R. B. Roncoli, "Lunar constants and models document," Jet
%    Propulsion Laboratory, California Institute of Technology, Tech.
%    Rep. JPL D-32296, 23 Sep. 2005. [Online]. Available: 
%    http://www.hq.nasa.gov/alsj/lunar cmd 2005 jpl d32296.pdf
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<3||isempty(timeCoordSys))
        timeCoordSys='TT';
    end
    
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

    if(~strcmp(timeCoordSys,'TDB')&&(nargin<4||isempty(deltaTTUT1)))
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
    
    EulerAngs=readJPLEphem(TDB1,TDB2,15,[]);
    R=Euler3Ang2RotMat(EulerAngs(1),EulerAngs(2),EulerAngs(3),'xyz');    
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
