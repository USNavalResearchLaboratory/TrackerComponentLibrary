function [year,month,day,hour,minute,second]=UTC2Cal(Jul1,Jul2,giveDayFrac)
%%UTC2Cal Convert universal coordinated times (UTC) given as a
%         two-part pseudo-Julian dates to a date in terms of the Gregorian
%         calendar in years, months, days, hours, minutes and seconds or
%         in terms of the Gregorian calendar in years, months, days and
%         the fraction of a day.
%
%INPUTS: Jul1, Jul2 Vectors or matrices of the two parts of Julian dates
%                   given in UTC. The units of the date are days. The full
%                   date is the sum of both terms. The date is broken into
%                   two parts to provide more bits of precision. It does
%                   not matter how the date is split.
%       giveDayFrac An optional boolean variable specifying whether the
%                   output should be as years months, days and a fraction
%                   of a day. If this parameter is omitted, the default is
%                   false. That means that the output will be given as
%                   years, months, days, hours, minutes, and seconds.
%
%OUTPUTS: Regardless of the value of giveDayFrac, the first three outputs
%         are the same. They are:
%           year A vector or matrix of years in the Gregorian calendar
%                under UTC time, one for for each Julian date.
%          month A vector or matrix of months in the Gregorian calendar
%                under UTC time, one for each Julian date. 1<=month<=12
%            day A vector or matrix of days in the Gregorian calendar
%                under UTC time, one for each Julian date. Days count from
%                1.
%          If giveDayFrac is omitted or is false, then three additional
%          outputs are
%           hour  A vector or matrix of hours in the Gregorian calendar
%                 under UTC time, one for each Julian date. 0<=hour<=23.
%          minute A vector or matrix of minutes in the Gregorian calendar
%                 under UTC time, one for each Julian date.
%                 0<=minute<=59.
%          second A vector or matrix of seconds in the Gregorian calendar
%                 under UTC time, one for each Julian date. This includes
%                 the possibility of a leap second on the day in question.
%          On the other hand, if giveDayFrac is true, then the one
%          additional output is
%          dayFrac A vector or matrix of values >=0 and <1 indicating the
%                  fraction of the day elapsed, one for each Julian date.
%
%The UTC date is only pseudo-Julian, because there is not a fixed number
%of seconds in a Julian day. The convention used in the IAU standard is
%that the Julian day matches the UTC day regardless of whether the UTC day
%is 86399, 86400 or 86401 SI seconds (depending on the presence of leap
%seconds).
%
%UTC began at 1960 January 1.0 (JD 2436934.5) and this function should not
%be called with an earlier date.
%
%The algorithm can be compiled for use in Matlab  using the
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[year,month,day,hour,minute,second]=UTC2Cal(Jul1,Jul2);
%or using
%[year,month,day,dayFrac]=UTC2Cal(Jul1,Jul2,true);
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

error('This function is only implemented as a mexed C or C++ function. Please run CompileCLibraries.m to compile the function for use.')

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
