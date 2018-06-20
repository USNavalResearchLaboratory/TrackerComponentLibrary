function [Jul1,Jul2]=Cal2UTC(year,month,day,hour,minute,second)
%%CAL2UTC Convert dates in terms of the Gregorian calendar in years,
%         months, days, hours, minutes and seconds with the time in
%         universal coordinated time (UTC) to a two-part pseudo-Julian
%         date in UTC.
%
%INPUTS: year A matrix of integer  years in the Gregorian calendar under
%             UTC time.
%       month A matrix of integer months in the Gregorian calendar under
%             UTC time. 1<=month<=12
%         day A matrix of integer days in the Gregorian calendar under UTC
%             time. Days count from 1.
%        hour A matrix of integer hours under the Gregorian calendar. UTC
%             time 0<=hour<=23
%      minute A matrix of integer minutes in the Gregorian calendar under
%             UTC time. 0<=minute<=59.
%      second A matrix of double floating point seconds in the Gregorian
%             calendar under UTC time. This is >=0 and normally less than
%             60, but can be a value less than 61 or 59 at the right hour
%             on a day with a leap second.
%
%Since the default format of numbers in Matlab is double-precision
%floating point numbers and one must go to lengths to use values that are
%differently formatted, it is assumed that all inputs, including the
%integer values, are double-precision floats.
%
%OUTPUTS: Jul1, Jul2 Matrices of the time as pseudo-Julian dates in UTC
%                    where each row/column corresponds to the values in
%                    the same row/column of the input matrices. Jul1
%                    corresponds to an integer number of days starting at
%                    midnight, which means that it ends in 0.5, since
%                    Julian days start at noon. Jul2 is the fraction of a
%                    day (midnight-to-midnight) after that.
%
%This is a mex wrapper for the function iauDtf2d in the International
%Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
%The UTC date is only pseudo-Julian, because there is not a fixed number
%of seconds in a Julian day. The convention used in the IAU standard is
%that the Julian day matches the UTC day regardless of whether the UTC day
%is 86399, 86400 or 86401 SI seconds (depending on the presence of leap
%seconds).
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[Jul1,Jul2]=Cal2UTC(year,month,day,hour,minute,second);
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
