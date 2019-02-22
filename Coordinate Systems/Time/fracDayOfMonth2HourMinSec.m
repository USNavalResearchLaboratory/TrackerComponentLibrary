function [hour,minute,second]=fracDayOfMonth2HourMinSec(year,month,day,dayFrac)
%%FRACDAYOFMONTH2HOURMINSEC Convert a fractional day of a month in the
%                   Gregorian calendar into the hour, minute, and (floating
%                   point) second of the day taking into account UTC
%                   leapseconds.
%
%INPUTS: year  A matrix of integer years (as Matlab doubles) in the
%              Gregorian calendar under UTC.
%       month  A matrix of integer months (as Matlab doubles) in the
%              Gregorian calendar under UTC (1-12).
%        day   A matrix of integer months (as Matlab doubles) in the
%              Gregorian calendar under UTC counting from 1.
%      dayFrac A matrix of fractional days in the Gregorian Calendar under
%              UTC (0<=dayFrac<1).
%
%OUTPUTS: hour A matrix of hours corresponding to time of day in the input
%              times (0<=hour<=23).
%       minute A matrix of minutes corresponding to time of day in the
%              input times (0<=minute<=59)
%       second A matrix of double floating point seconds in the Gregorian
%              calendar under UTC time corresponding to the input day
%              fraction. This is >=0 and normally less than 60, but can
%              be a value less than 61 or 59 at the right hour on a day
%              with a leap second.
%
%The function leapSecsOnDay is used to determine whether a day contains
%leap seconds. After that, the rest is fairely straightforward.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numValues=length(year(:));

%Allocate the return values.
hour=zeros(size(year));
minute=zeros(size(year));
second=zeros(size(year));

for curVal=1:numValues
    if(dayFrac(curVal)>=1)
       warning('The time is after the end of the day.')
    end
    
    if(dayFrac(curVal)<0)
        error('Invalid day fraction entered')
    end
    
    secInDay=60*60*24+leapSecsOnDay(year(curVal),month(curVal),day(curVal));
    numSec=dayFrac(curVal)*secInDay;

    %The min takes care of the leapsecond case where the second number can hit
    %60 in the final hour of the day in UTC time.
    hour(curVal)=min(fix(numSec/(60*60)),23);

    %Minutes are found with residual seconds after getting rid of seconds in
    %hours, which will never hold one of the leapseconds.
    numSec=numSec-60*60*hour(curVal);
    %The min takes care of the leapsecond case in the final hour
    minute(curVal)=min(fix(numSec/(60)),59);

    %Residual seconds, including the leapsecond.
    second(curVal)=numSec-minute(curVal)*60;
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
