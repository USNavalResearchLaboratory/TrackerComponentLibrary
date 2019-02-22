function [UTC1,UTC2]=UT12UTC(UT11,UT12,deltaUTCUT1)
%%UT12UTC Convert from UT1, which is a nonuniform timescale based on the
%         rotation rate of the Earth, to Coordinated Universal Time (UTC),
%         which is the "standard" time that pretty much every country uses
%         (the time zone is that of the United Kingdom) including
%         leapseconds. Due to the use of leapseconds, UT1 is always within
%         0.9 seconds of UTC.
%
%INPUTS: Jul1,Jul2 Two parts of a Julian date given in UTC. The units of
%                  the date are days. The full date is the sum of both
%                  terms. The date is broken into two parts to provide
%                  more bits of precision. It does not matter how the date
%                  is split.
%      deltaUTCUT1 An optional parameter specifying the offset between UTC
%                  and UT1 in seconds. If this parameter is omitted, then
%                  the value of the function getEOP will be used and
%                  iterations will be performed, because the correct
%                  tabulated value is given in UTC, not UT1.
%
%OUTPUTS: Jul1,Jul2 Two parts of a pseudo-Julian date in UTC with units of
%                   Julian days.
%
%This essentially just adds deltaUTCUT1 from the given UTC date, or
%looks up deltaUTCUT1 using the getEOP function and then adds
%deltaUTCUT1. The "looking up" step involves an interation, since the value
%of deltaUTCUT1 is tabulated by UTC, which is unknown. To know the peroper
%amount to subtract from the pseudo-Julian day, the function leapSecsOnDay
%is used to determine whether the day contains a leapsecond. As that
%function is also tabulated by UTC, an iteration is necessary even if
%deltaUTCUT1 is given. Unlike the function iauUtcut1 in the International
%Astronomical Union's Standards of Fundamental Astronomy library, this does
%not try to deal with deltaUTCUT1 being off by a second near a leapsecond
%as it was noticed that the "corrections" made in the library make the pair
%of UTC2UT1 and UT12UTC inconsistent by a second when it is near a
%leapsecond.
%
%Many temporal coordinate systems standards are compared in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report, Naval
%    Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016, 173
%    pages.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The initial estimate of UTC is just UT1. This will be within 0.9 seconds
%of the correct value.
UTC1=UT11;
UTC2=UT12;

for curIter=1:2
    if(nargin<3)
        [~,~,deltaUTCUT1]=getEOP(UTC1,UTC2);
    end
    
    %Find the number of seconds in the day.
    [year,month,day]=UTC2Cal(UTC1,UTC2,true);
    numSecInDay=24*60*(60+leapSecsOnDay(year,month,day));

    %Add preserving precision.
    if(UTC1<UTC2)
        UTC2=UT11+deltaUTCUT1/numSecInDay;
        UTC1=UT12;
    else
        UTC2=UT12+deltaUTCUT1/numSecInDay;
        UTC1=UT11;
    end

    %Put the fractional part into UTC2 and the integer part into UTC1.
    frac1=UTC1-fix(UTC1);
    frac2=UTC2-fix(UTC2);
    sumVal=frac1+frac2;
    sumInt=fix(sumVal);
    sumFrac=sumVal-sumInt;
    UTC1=fix(UTC2)+fix(UTC1)+sumInt;
    UTC2=sumFrac;
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
