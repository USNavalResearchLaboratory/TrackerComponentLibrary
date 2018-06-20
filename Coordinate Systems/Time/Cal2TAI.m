function [Jul1,Jul2]=Cal2TAI(year,month,day,hour,minute,second)
%%CAL2TT  Convert dates in terms of the Gregorian calendar in years,
%          months, days, hours, minutes and seconds with the time in
%          universal coordinated time (UTC) to a two-part Julian
%          date in international atomic time (TAI).
%
%INPUTS:   year    A matrix of integer  years in the Gregorian
%                  calendar under UTC time.
%          month   A matrix of integer months in the Gregorian calendar
%                  under UTC time. 1<=month<=12
%          day     A matrix of integer days in the Gregorian calendar
%                  under UTC time. Days count from 1.
%          hour    A matrix of integer hour under the Gregorian calendar.
%                  UTC time 0<=hour<=23
%          minute  A matrix of integer minutes in the Gregorian calendar
%                  under UTC time. 0<=minute<=59.
%          second  A matrix of double floating point second in the
%                  Gregorian calendar under UTC time. This is normally
%                  less than 60, but can be a value less than 61 or 59 at
%                  the right hour on a day with a leap second.
%
%OUTPUTS:   Jul1, Jul2  Matrices of the time as pseudo-Julian dates in TAI
%                       where each row/column corresponds to the values in
%                       the same row/column of the input matrices.
%
%This function just calls the function Cal2UTC and then UTC2TAI.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[Jul1,Jul2]=Cal2UTC(year,month,day,hour,minute,second);
[Jul1,Jul2]=UTC2TAI(Jul1,Jul2);
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
