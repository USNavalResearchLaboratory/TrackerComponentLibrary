function numLS=leapSecsOnDay(year,month,day)
%%LEAPSECSONDAY Returns the number of leap seconds in a day for dates from
%               1972 onward. Prior to 1972, non-integer leap second numbers
%               were used, meaning that the time of day matters. From 1972
%               onward, one could have a positive or negative leap second,
%               meaning that the last minute of the day in coordinated
%               universal time (UTC) might contain 59 or 61 seconds.
%
%INPUTS: year    An integer  year in the Gregorian calendar under UTC
%                time. year>= 1960 when UTC started.
%       month    An integer month in the Gregorian calendar under UTC
%                time. 1<=month<=12
%        day     An integer day in the Gregorian calendar under UTC
%                time. Days count from 1.
%
%OUTPUTS: numLS  The number of leap seconds on the given day. This can be
%                0, -1, or 1.
%
%This function calls the cumLeapSec function for the day in question and
%for one day later to see what the difference is. The function cumLeapSec
%relies on iaudat in the International Astronomical Union's (IAU)
%Standard's of Fundamental Astronomy library. Note that leap second
%information is not available for future dates as they can not be
%accurately predicted.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(year<1972)
        error('Dates before 1972 have a more complicated leap second methodology than is used here');
    end

%Determine the number of leap on the date.
    numSec=cumLeapSec(year,month,day);
    [Jul1,Jul2]=Cal2UTC(year,month,day,0,0,0);
    %Advance the UTC day count by one day.
    Jul1=Jul1+1;
    [year,month,day]=UTC2Cal(Jul1,Jul2,true);
    %Determine the number of leap seconds on day later.
    numSecNextDay=cumLeapSec(year,month,day);
    
    numLS=numSecNextDay-numSec;
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
