function [rGCRS,rITRS]=solarBodyVec(Jul1,Jul2,timeCoordSys,objectNumber,obsState,observerCoordSys,deltaTTUT1,xpyp,dXdY)
%%SOLARBODYVEC Get a vector from a near-Earth object to a given solar body,
%              such as a planet or the Sun, including basic light-time and
%              aberration corrections. Refraction corrections must be
%              applied separately. The observer can be in the ITRS (WGS-84)
%              or GCRS coordinate systems.
%
%INPUTS: Jul1,Jul2 Two parts of a Julian date given in the time system
%              specified by the second argument. The units of the date are
%              days. The full date is the sum of both terms. The date is
%              broken into two parts to provide more bits of precision. It
%              does not matter how the date is split.
% timeCoordSys A parameter specifying the coordinate system of the Julian
%              date given by Jul1,Jul2. Possible values are
%              'UTC' Universal coordinated time.
%              'TT'  Terrestrial time.
%              'TDB' Barycentric dynamical time (can only be used if
%                    algorithm =0 or 1)
%              'TAI' International atomic time.
%              'TCG' Geocentric coordinate time.
%              'GPS' GPS time.
%              If the input is not given in TDB, then it must be converted
%              to TDB. If the input is given in TDB but the observer's
%              location is in the ITRS, or rITRS is desired on the output,
%              then it must be converted to TT. In computing TDB or TT,
%              the observer's clock is assumed to be collocated with the
%              center of the Earth --thus, no topocentric corrections are
%              applied. This might result in an often insignificant loss of
%              precision. Also, the functions TDB2TT and TT2TDB are used,
%              which are less accurate than using ephemerides to determine
%              the offset between TDB and TT.
% objectNumber An integer indicating which celestial object is being
%              observed. Possible values are:
%              1  Mercury
%              2  Venus 
%              3  Earth (geocenter)
%              4  Mars (system barycenter)
%              5  Jupiter (system barycenter)
%              6  Saturn (system barycenter)
%              7  Uranus (system barycenter)
%              8  Neptune (system barycenter)
%              9  Pluto (system barycenter)
%              10 Moon
%              11 Sun
%              12 Solar System Barycenter
%              13 Earth-Moon Barycenter
%     obsState The 6X1 location and velocity vector an observer in
%              Cartesian coordinates (in meters) in the coordinate system
%              selected by the next argument. The format is
%              [x;y;z;xdot;ydot;zdot].
% stateCoordSys The coordinate system that obsState is in. The
%              possibilities are:
%              'GCRS' The geocentric celestial reference system.
%              'ITRS' The international terrestrial reference system.
%   deltaTTUT1 An optional parameter specifying the difference between TT
%              and UT1 in seconds. Values of this parameter can be obtained
%              from
%              http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%              or 
%              http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%              If this parameter is omitted or an empty matrix is passed,
%              then the value provided by the function getEOP will be used
%              instead.
%         xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
%              including the effects of tides and librations. If this
%              parameter is omitted or if an empty matrix is passed, the
%              value from the function getEOP will be used.
%         dXdY dXdY=[dX;dY] are the celestial pole offsets with respect to
%              the IAU 2006/2000A precession/nutation model in radians If
%              this parameter is omitted or an empty matrix is passed, the
%              value from the function getEOP will be used.
%
%OUTPUTS: rICRS A 3X1 vector in Cartesian GCRS coordinates in meters from
%               the observer to the object, accounting for light time and
%               aberration corrections.
%         rITRS A vector in Cartesian ITRS coordinates in meters from the
%               observer to the object with any light-time and aberration
%               corrections as requested. 
%
%The relative location of the celestial body is given by the readJPLEphem
%function, which reads the National Aeronautics and Space Administration
%(NASA) Jet Propulsion Laboratory's (JPL's) digital ephemerides. Light-time
%and aberration corrections are applied. Those effects are defined in
%Chapter 7 of [1]. Note that unlike the function iauAb in the IAU's
%Standards of Fundamental Astronomy library, no attempt is made to correct
%for the Sun's gravitational potential, which adds an error of up to 4
%microarcseconds. Additionally, gravitational deflection from the Sun and
%other planets is not taken into account.
%
%If getEOP has to be used, deltaTTUT1 is not provided and the time is given
%in TDB, then the terrestrial time (TT) needed for the lookup is
%approximated using TDB.
%
%REFERENCES:
%[1] S. E. Urban and K. P.Seidelmann, Eds.,Explanatory Supplement to the
%    Astronomical Almanac, 3rd ed. Mill Valley, CA: University Science
%    Books, 2013.
%
%August 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<9)
    dXdY=[];
end

if(nargin<8)
    xpyp=[];
end

if(nargin<7)
    deltaTTUT1=[];
end

%The time needs to be in TDB. Perform a conversion, if it does not depend
%on the chosen coordinate system. If given TT or UTC, convert to TT and
%wait until the coordinate system for the location of the observer is
%known.
switch(timeCoordSys)
    case 'GPS'
        [TT1,TT2]=GPS2TT(Jul1,Jul2);
    case 'TCG'
        [TT1,TT2]=TCG2TT(Jul1,Jul2);
    case 'TT'
    case 'TDB'
        %If the time is given in TDB, we will only need TT if the observer
        %state is given in the ITRS.
        TDB1=Jul1;
        TDB2=Jul2;
    case 'TAI'
        [TT1,TT2]=TAI2TT(Jul1,Jul2);
    case 'UTC'
        [TT1,TT2]=UTC2TT(Jul1,Jul2);
    otherwise
        error('An unsupported coordinate system for the time is given.')
end

if(isempty(dXdY)||isempty(xpyp)||isempty(deltaTTUT1))
    if(~exist('TT1','var'))
        if(~isempty(deltaTTUT1))
            [TT1t,TT2t]=TDB2TT(TDB1,TDB2,deltaTTUT1);
        else
            %If TDB is given and deltaTTUT1 is not given, the approximate
            %TT with TDB when querying the EOPs.
            TT1t=TDB1;
            TT2t=TDB2;
        end
    else
        TT1t=TT1;
        TT2t=TT2;
    end

    [UTC1,UTC2]=TT2UTC(TT1t,TT2t);
    [xpyp1,dXdY1,~,deltaTTUT11]=getEOP(UTC1,UTC2);
    
    if(isempty(xpyp))
       xpyp=xpyp1; 
    end
    
    if(isempty(dXdY))
       dXdY=dXdY1; 
    end
    
    if(isempty(deltaTTUT1))
       deltaTTUT1=deltaTTUT11; 
    end
end

if(~strcmp(timeCoordSys,'TDB'))
    [TDB1,TDB2]=TT2TDB(TT1,TT2,deltaTTUT1);
end

switch(observerCoordSys)
    case 'GCRS'%No conversion of the state is needed.
        obsPosGeo=obsState(1:3);
        obsVelGeo=obsState(4:6);
    case 'ITRS'
        %Convert the position and velocity components of obsState into
        %GCRS coordinates. The includes adding in the velocity due to the
        %Earth's rotation rate. This requires that we have terrestrial
        %time.
        if(strcmp(timeCoordSys,'TDB'))
            [TT1,TT2]=TDB2TT(TDB1,TDB2,deltaTTUT1);
        end
        
        obsGCRS=ITRS2GCRS(obsState,TT1,TT2,deltaTTUT1,xpyp,dXdY);

        obsPosGeo=obsGCRS(1:3);
        obsVelGeo=obsGCRS(4:6);
    otherwise
       error('Unsupported Observer coordinate system chosen.')
end

%Perform three iterations of the light-time correction (ignoring
%gravitational deflection of the light due to other planets or the Sun).
c=Constants.speedOfLight;%Units of meters per second.
LTEst=0;
%Assume that it will converge in three iterations.
for curIter=1:3
    if(TDB1>TDB2)
       TDB1Cur=TDB1;
       TDB2Cur=TDB2-LTEst;
    else
       TDB1Cur=TDB1-LTEst;
       TDB2Cur=TDB2;
    end
    
    state=readJPLEphem(TDB1Cur,TDB2Cur,objectNumber,3);
                       
%The returned state is a vector of 3D position and velocity with units
%of meters and meters per second. Compute the time it would take light to
%travel from the object to the observer
    relPos=state(1:3)-obsPosGeo;
    LTEst=norm(relPos)/c;%The units are seconds.
    %Convert to Julian days.
    LTEst=LTEst/86400;
end

%Next, apply a special relativistic correction to the postion of the
%object. For this, we need the instantaneous velocity vector of the
%state coordinate system in the BCRS coordinate system when the measurement
%is taken.
EarthStateBCRS=readJPLEphem(TDB1,TDB2,3,12);
        
%The velocity of the observer in the BCRS coordinate system is the velocity
%of the Earth plus the velocity of the observer with respect to the Earth.
%The addition is performed in a special relativistically correct manner.
obsVelBCRS=relVecAdd(EarthStateBCRS(4:6),obsVelGeo);

%Get the distance from the origin of the state coordinate system to
%the Sun.
SunToEarthState=readJPLEphem(TDB1,TDB2,3,11);

%The distance from the observer to the Sun.
sunDist=norm(SunToEarthState(1:3)-obsPosGeo);

%Perform the aberration correction.
rGCRS=aberCorr(relPos,obsVelBCRS,sunDist);

%If the position in the ITRS is desired.
if(nargout>1)
    if(~exist('TT1','var'))
        [TT1,TT2]=TDB2TT(TDB1,TDB2,deltaTTUT1);
    end
    rITRS=GCRS2ITRS(rGCRS,TT1,TT2,deltaTTUT1,xpyp,dXdY);
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
