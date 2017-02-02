function [yearFrac,secOfYear]=TT2YearFrac(Jul1,Jul2)
%%TT2YEARFRAC A function for converting a terrestial time (TT) given as a
%              two-part pseudo-Julian to a decimal year in the Gregorian
%              calendar in universal coordinated time (UTC), accounting for
%              leap years and leap seconds. The fractional part of the year
%              is the ratio of the number of seconds elapsed in the year to
%              the total number of seconds in the year. This function also
%              returns the total number of seconds in the year. This is
%              only valid from 1960 onward, as UTC was established in 1960.
%
%INPUTS: Jul1, Jul2  Vectors or matrices of the two parts of Julian dates
%                    given in TT. The units of the date are days. The full
%                    date is the sum of both terms. The date is broken into
%                    two parts to provide more bits of precision. It does
%                    not matter how the date is split.
%
%OUTPUTS: yearFrac   A matrix of the year as a fractional date
%                    corresponding to the input dates. For example, 2005.25
%                    means that the year is 2005 and 1/4 of the seconds in
%                    the year have passed.
%        secOfYear   A matrix of the fractional seconds from the start of
%                    the year as implied by the given dates.
%
%This function just calls the function TT2UTC followed by the function
%UTC2YearFrac.
%
%Note that the type of fractional years returned by this function should
%not be confused with Julian epochs, which are very similar and can be
%obtained using the function JulDate2JulEpoch.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[Jul1,Jul2]=TT2UTC(Jul1,Jul2);
[yearFrac,secOfYear]=UTC2YearFrac(Jul1,Jul2);
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
