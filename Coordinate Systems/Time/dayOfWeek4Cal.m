function dayOfWeek=dayOfWeek4Cal(year,month,day)
%%DAYOFWEEK4CAL Determine the day of the week given a date in terms of the
%               Gregorian calendar in years, months, and days in universal
%               coordinated time (UTC).
%
%INPUTS:   year    A matrix of integer  years in the Gregorian
%                  calendar under UTC time.
%          month   A matrix of integer months in the Gregorian calendar
%                  under UTC time. 1<=month<=12
%          day     A matrix of integer days in the Gregorian calendar
%                  under UTC time. Days count from 1.
%
%OUTPUTS: dayOfWeek A matrix of integers indicating the day of the week for
%                   each date given. 0=Sunday, 1=Monday, 2=Tuesday,
%                   3=Wednesday, 4=Thursday, 5=Friday,6=Saturday.
%
%Since the days repeat in a 7 day cycle, this just turns the date into a
%pseudo-Julian date in UTC, which is in terms of days, adds 5, truncates it
%(because Julian dates start at noon and weekdays start at midnight) and
%takes the modulo 7. 
%
%August 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    zeroMat=zeros(size(year));
    [Jul1,Jul2]=Cal2UTC(year,month,day,zeroMat,zeroMat,zeroMat);
    %The fix function is because the days begin at noon in pseudo Julin
    %dates, but weekdays begin at midnight.
    dayOfWeek=mod(fix(Jul1+Jul2-5),7);
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
