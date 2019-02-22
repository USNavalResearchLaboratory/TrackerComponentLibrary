function [timeZoneInfo,displayName]=getTimeZoneInfo(timeZoneName,Jul1,Jul2)
%GETTIMEZONEINFO Get information about a time zone. If a date is provided,
%                one can obtain an offset from UTC taking daylight savings
%                time into account, if relevant. Note, however, that many
%                countries only add a leap second at their local midnight
%                and not when UTC does, so the offsets might be off by a
%                second around when a leap second is added.
%
%INPUTS: timeZoneName The name of the time zone. This is one of the short
%                   names returned by getTimeZoneList(). A meaningless
%                   input will cause UTC to be returned.
%           Jul1,Jul2 Optionally, the time as a pseudo-Julian date in UTC.
%                   The units of the date are days. The full date is the
%                   sum of both terms. The date is broken into two parts to
%                   provide more bits of precision. It does not matter how
%                   the date is split. If this input is provided, then
%                   daylight savings time information can be provided on
%                   the output.
%
%OUTPUTS: timeZoneInfo A 2X1 or 4X1 (if a UTC date is provided) providing
%                   information on the time zone. timeZoneInfo(1) is the
%                   offset from UTC in seconds for time time zone, without
%                   daylight savings time. timeZoneInfo(2) is true if the
%                   time zone is listed as having ever supported daylight
%                   savings time or if it is listed as ever supporting it
%                   in the future. timeZoneInfo(3) is the offset in seconds
%                   from UTC at the given date and timeZoneInfo(4) is true
%                   if the time zone is in daylight savings time at the
%                   given date.
%       displayName A more detailed version of timeZoneName.
%
%This function just calls the appropriate Java commands to obtain time zone
%information.
%
%EXAMPLE:
%For the Easten Time Zone of the United States on 11 January 2016, one can
%use
% [UTC1,UTC2]=Cal2UTC(2016,1,11,0,0,0);
% [timeZoneInfo,displayName]=getTimeZoneInfo('EST5EDT',UTC1,UTC2)
%One will see that daylight savings time is honored, but is not active at
%the given time, so both offsets are -18000 second (=-5 hours).
%
%May 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

timeZone=java.util.TimeZone.getTimeZone(timeZoneName);

displayName=char(timeZone.getDisplayName());

if(nargin>1&&~isempty(Jul1))
    timeZoneInfo=zeros(4,1);
    timeZoneInfo(1)=timeZone.getRawOffset()/1e3;
    timeZoneInfo(2)=timeZone.observesDaylightTime() || timeZone.useDaylightTime();
    
    [TT1,TT2]=UTC2TT(Jul1,Jul2);
    [UTCRef1,UTCRef2]=Cal2UTC(1970,1,1,0,0,0);
    [TTRef1,TTRef2]=UTC2TT(UTCRef1,UTCRef2);

    %It should be that TT1>TT2 and TTRef1>TTRef2 due to how the conversion
    %functions work.
    %There are 86400 seconds in all Julian days in terrestrial time. The
    %multiplication by 1000 turns it into milliseconds.
    timeDiffMilliSeconds=int64((TT1-TTRef1)*86400*1000)+int64((TT2-TTRef2)*86400*1000);

    timeZoneInfo(3)=timeZone.getOffset(timeDiffMilliSeconds)/1e3;
    timeZoneInfo(4)=timeZoneInfo(1)~=timeZoneInfo(3);%Must currently be in DST.
else
    timeZoneInfo=zeros(2,1);
    timeZoneInfo(1)=timeZone.getRawOffset()/1e3;
    timeZoneInfo(2)=timeZone.observesDaylightTime() || timeZone.useDaylightTime();
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
