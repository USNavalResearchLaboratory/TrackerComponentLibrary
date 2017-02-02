function [month,dayOfMonth]=dayOfYear2MonthDay(year,dayOfYear)
%%DAYOFYEAR2MONTHDAY Given a count of the day of the year (starting from 1)
%                    find the corresponding month of the year and the day
%                    of that month, accounting for leap years. Everything
%                    is assumed to be under the Gregorian calendar.
%
%INPUTS: year       A matrix of integer years in the Gregorian calendar.
%        dayOfYear  A matrix of integer days of the year in the Gregorian
%                   calendar starting from 1. 
%                   
%OUTPUTS: month A matrix of the corresponding months of the year for each
%               year-day pair inputted (1-12).
%    dayOfMonth A matrix of integer days of the month corresponding to each
%               day of the year value.
%
%This function basically just uses the isLeapYear function to determine if
%it is a leap year and then counts the number of cumulative days in each
%month, seeing in which month the dayOfYear falls.
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(any(dayOfYear<=0|dayOfYear~=fix(dayOfYear)))
        error('dayOfYear must be a positive integer.')
    end

    lastDaysOfMonthLeap=cumsum([31;29;31;30;31;30;31;31;30;31;30;31]);
    lastDaysOfMonthNonLeap=cumsum([31;28;31;30;31;30;31;31;30;31;30;31]);

    %Allocate space for return values.
    numEntries=length(year(:));
    month=zeros(size(year));
    dayOfMonth=zeros(size(year));  
    for curEntry=1:numEntries
        if(isLeapYear(year(curEntry)))
            lastDaysOfMonth=lastDaysOfMonthLeap;
        else
            lastDaysOfMonth=lastDaysOfMonthNonLeap;
        end

        month(curEntry)=find(dayOfYear(curEntry)<=lastDaysOfMonth,1);
        if(isempty(month(curEntry)))
            error('Invalid day of year provided')
        end

        if(month(curEntry)>1)
            dayOfMonth(curEntry)=dayOfYear(curEntry)-lastDaysOfMonth(month(curEntry)-1);
        else
            dayOfMonth(curEntry)=dayOfYear(curEntry);
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
