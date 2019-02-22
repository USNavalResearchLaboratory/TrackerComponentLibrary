function numSec=cumLeapSec(year,month,day,dayFrac)
%%CUMLEAPSEC Get the cumulative number of leap seconds up to a given date
%            for one or more dates given in matrices. The result is the
%            difference between atomic international time (TAI) and
%            universal coordinated time (UTC). The result is not always an
%            integer prior to 1972, because non whole second adjustments
%            were added back then. The results are only valid from the
%            start of UTC, which was on 1960 January 1.0.
%
%INPUTS: year A matrix of integer years in the Gregorian calendar under
%             UTC time.
%       month A matrix of integer months in the Gregorian calendar under
%             UTC time. 1<=month<=12
%         day A matrix of integer days in the Gregorian calendar under
%             UTC time. Days count from 1.   
%     dayFrac A parameter that is only neded for years prior to 1972,
%             because of the unusual manner in which time corrections were
%             applied back then. The fraction of a day must be between 0
%             and 1 for each element in the matrix. If this parameter is
%             omitted, then a value of zero is used for all dates.
%
%This function is just a wrapper for the function iaudat in the
%International Astronomical Union's (IAU) Standard's of Fundamental
%Astronomy library. If the year in question is too far into the future
%(from the version of the IAU's function) the function will issue a
%warning message as it does not hold a record of future leap seconds and
%leap seconds are not accurately precitable.
%
%The algorithm can be compiled for use in Matlab  using the
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%numSec=cumLeapSec(year,month,day);
%for 1972 onwards and the format
%numSec=cumLeapSec(year,month,day,dayFrac);
%for dates prior to 1972.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
