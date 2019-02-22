function val=isLeapYear(year)
%%ISLEAPYEAR Returns true if the year in question is a leap year (has 366
%            days) and 1 if the year has only 365 days.
%
%INPUTS: year A matrix of integer year numbers in the Gregorian calendar.
%             These should be the year 1753 or later.
%
%OUTPUTS: val A matrix of boolean values corresponding to the years in year
%             that are true if the year is a leap year and false if the
%             year is not a leap year.
%
%The Gregorian calendar was introduced by Pope Gregory XIII on the 24th of
%February, 1582 to replace the Julian calendar. To deal with the loss of
%days due to the inaccuracy of the Julian calendar, ten days were deleted
%from October of that year, so October 4th was followed by October 15th.
%
%However, the calendar was rejected by most protestant countries for many
%years. It was only in 1752 that England accepted the Gregorian calendar,
%dopping 11 days from September. Due to the ambiguity of dates (at last
%when considering European and American history) prior to 1753, this
%function issues an error when one tries to determine whether a year prior
%to 1753 was a leap year.
%
%More information on the gregorian calendar and leap years is given in [1].
%
%February has 29 days if the year is divisible by 4, unless the year is
%divisible by 100 unless the year is divisible by 400.
%
%REFERENCES:
%[1] E. G. Richards, "Calendars," in Explanatory Supplement to the
%    Astronomical Almanac, 3rd ed., S. E. Urban and P. K. Seidelmann, Eds.
%    Mill Valley, CA: University Science Books, 2013, ch. 15.
%
%January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    %Allocate soace.
    val=zeros(size(year));
    numYear=length(year(:));
    
    for curYear=1:numYear
        if(year(curYear)<1753)
           error('Years prior to 1753 are not supported') 
        end

        if(mod(year(curYear),4)==0)
            if(mod(year(curYear),100)==0&&mod(year(curYear),400)~=0)
                val(curYear)=false;%Not a leap year.
            else
                val(curYear)=true;%It is a leap year.
            end
        else%Not a leap year.
            val(curYear)=false;
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
