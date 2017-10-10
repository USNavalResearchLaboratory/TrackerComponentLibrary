function [LATRad,equationOfTime,LAT1,LAT2]=TT2LAT(TT1,TT2,rObsITRS,deltaTTUT1,xpyp,dXdY)
%%TT2LAT Convert from terrestrial time to local apparent solar time (LAT).
%        This is the time that one would read on a sundial that is
%        stationary in the Terrestrial intermediate Reference System (TIRS)
%        taking into account light-time and aberration but not atmospheric
%        refraction. Unlike the International Terrestrial Reference System
%        (ITRS), the TIRS does not take into account polar motion; that is,
%        the motion of the rotation axis of the Earth with respect to the
%        crust of the Earth. Since the polar motion makes the rotation axis
%        wobble with respect to the crust, an observer at a constant
%        location in the ITRS will move over time as seen in the TIRS. The
%        location of the Sun is determined via the function readJPLEphem
%        whereby TDB is approxmated using TT2TDB. The equation of time is
%        also returned.
%
%INPUTS: TT1,TT2 Two parts of a Julian date given in TT. The units of the
%                date are days. The full date is the sum of both terms. The
%                date is broken into two parts to provide more bits of
%                precision. It does not matter how the date is split.
%       rObsITRS The 3X1 location of the observer in the International
%                Terrestrial Reference System (ITRS). Only the direction
%                matters, not the magnitude.
%     deltaTTUT1 An optional parameter specifying the difference between TT
%                and UT1 in seconds. This information can be obtained from
%                http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%                or 
%                http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%                If this parameter is omitted or if an empty matrix is
%                passed, then the value provided by the function getEOP
%                will be used instead.
%           xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
%                including the effects of tides and librations. If this
%                parameter is omitted or if an empty matrix is passed,
%                the value from the function getEOP will be used.
%           dXdY dXdY=[dX;dY] are the celestial pole offsets with respect
%                to the IAU 2006/2000A precession/nutation model in
%                radians. If this parameter is omitted, the value from the
%                function getEOP will be used.
%
%OUTPUTS:  LATRad The local apparent solar time (LAT) in radians. 0 radians
%                 is solar midnight. pi radians is solar noon. The mapping
%                 is 2*pi radians for 24 hours. Thus, LATRad*(24/(2*pi))
%                 gives the time of day in hours.
%  equationOfTime The equation of time in radians. This is the difference
%                 between LAT and local mean solar time (LMT).
%       LAT1,LAT2 Two parts of the local apparent solar time in Julian
%                 days. The date is split so that LAT1 is the integer part
%                 and LAT2 is the fractional part. Like UT1, a zero
%                 fractional part of the day corresponds to noon, not
%                 midnight.
%
%The local apparent solar time (LAT) is the local hour angle of the Sun
%(given as time) plus 12 hours. LAT is discussed in Chapter 1.2 of [1]. The
%"apparent" in the name means that light-time and aberration were taken
%into account (The glossary of the book defined "apparent". The coordinate
%system used for the computations is theTIRS, not the ITRS.
%
%REFERENCES:
%[1] S. E. Urban and K. P.Seidelmann, Eds.,Explanatory Supplement to the
%    Astronomical Almanac, 3rd ed. Mill Valley, CA: University Science
%    Books, 2013.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Get any Earth orientation parameters that were not provided.
if(nargin<6||isempty(deltaTTUT1)||isempty(xpyp)||isempty(dXdY))
    [UTC1,UTC2]=TT2UTC(TT1,TT2);
    [xpypNew,dXdYNEW,~,deltaTTUT1NEW]=getEOP(UTC1,UTC2);
    
    if(nargin<4||isempty(deltaTTUT1))
        deltaTTUT1=deltaTTUT1NEW;
    end
    
    if(nargin<5||isempty(xpyp))
        xpyp=xpypNew;
    end
    
    if(nargin<6||isempty(dXdY))
        dXdY=dXdYNEW;
    end
end

rOBSTIRS=ITRS2TIRS(rObsITRS,TT1,TT2,xpyp);
rObsSpher=Cart2Sphere(rOBSTIRS);
%Get the location of the Sun including light-time and aberration
%corrections.

%Sun Position with respect to Earth.
[TDB1,TDB2]=TT2TDB(TT1,TT2);
SunGCRSPosVel=readJPLEphem(TDB1,TDB2,11,3);
rSunTIRS=GCRS2TIRS(SunGCRSPosVel(1:3),TT1,TT2,deltaTTUT1,xpyp,dXdY);

rSunSpher=Cart2Sphere(rSunTIRS);

%The hour angle of the Sun.
HASun=-rSunSpher(2);

%The LOCAL hour angle of the Sun.
LHASun=HASun+rObsSpher(2);

%The local apparent solar time is the local hour angle of the sun plus 12
%hours given as an angle (24 hours is 2*pi radians). The 12 hours is
%because 0 for LAT is midnight, not noon as is the case for UT1.
LATRad=wrapRange(LHASun+pi,-pi,pi);

%If a full Julian date is desired, then find the equation of time. This is
%the difference between LMT and LAT. This is going to come from the
%fractional part of the day. Given that difference, it can be added to the
%Julian date for LMT to get the LAT Julian date.
if(nargout>1)
    %Get local mean time.
    [LMTRad,LMT1,LMT2]=TT2LMT(TT1,TT2,rObsITRS,deltaTTUT1,xpyp);

    equationOfTime=wrapRange(LATRad-LMTRad,-pi,pi);

    %LMT2 will be smaller than LMT1, so the equation of time is added to
    %that.
    LAT1=LMT1;
    %The (1/(2*pi)) factor is because there are exactly 2*pi radians per
    %day and LAT2 is denominated in Julian days.
    LAT2=LMT2+equationOfTime*(1/(2*pi));

    %Put the fractional part into LMT2 and the integer part into LMT1.
    frac1=LAT1-fix(LAT1);
    frac2=LAT2-fix(LAT2);
    sumVal=frac1+frac2;
    sumInt=fix(sumVal);
    sumFrac=sumVal-sumInt;
    LAT1=fix(LAT2)+fix(LAT1)+sumInt;
    LAT2=sumFrac;
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
