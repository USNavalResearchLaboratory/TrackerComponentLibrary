function numSecs=numSecInYear(year)
%%NUMSECINYEAR Find the number of seconds in a particular Gregorian
%              calendar year in universal coordinated time (UTC) from 1960
%              onwards, accounting for leap years and leap seconds.
%
%INPUTS: year A matrix  of integer years in the Gregorian calendar. The
%             results are only valid from the start of UTC, which was at
%             the start of 1960.
%
%OUTPUTS: numSecs A matrix of the number of seconds in the given years.
%
%The number of leap seconds in a year is determined by taking the
%difference of the cumulative number of leap seconds at the star and end of
%the year using the cumLeapSec function. The number of days in the year is
%then determined with the isLeapYear function. 24 hours per day,
%60 minutes per hour, 60 seconds per minute, times the number of days in
%the year is the number of seconds per year without leap seconds to which
%the leap seconds are added.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

leapSecs=cumLeapSec(year+1,ones(size(year)),ones(size(year)))-cumLeapSec(year,ones(size(year)),ones(size(year)));
leapYears=isLeapYear(year);

numYear=length(year(:));
numSecs=zeros(size(year));

for curYear=1:numYear
    if(leapYears(curYear)==true)
        numSecs(curYear)=leapSecs(curYear)+24*60*60*366;
    else
        numSecs(curYear)=leapSecs(curYear)+24*60*60*365;
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
