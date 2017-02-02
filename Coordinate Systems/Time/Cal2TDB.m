function [Jul1,Jul2]=Cal2TDB(year,month,day,hour,minute,second,deltaT,clockLoc)
%%CAL2TDB  Convert dates in terms of the Gregorian calendar in years,
%          months, days, hours, minutes and seconds with the time in
%          universal coordinated time (UTC) to a two-part Julian
%          date in barycentric dynamic time (TDB).
%
%INPUTS:   year    A year in the Gregorian calendar under UTC time.
%          month   A month in the Gregorian calendar under UTC time.
%                  1<=month<=12
%          day     An integer day in the Gregorian calendar under UTC time.
%                  Days count from 1.
%          hour    An integer hour under the Gregorian calendar.
%                  UTC time 0<=hour<=23
%          minute  Aninteger minutes in the Gregorian calendar
%                  under UTC time. 0<=minute<=59.
%          second  A double floating point second in the Gregorian calendar
%                  under UTC time. This is normally less than 60, but can
%                  be a value less than 61 or 59 at the right hour on a
%                  day with a leap second.
%        deltaT    An optional parameter specifying the offset between TT
%                  and UT1 in seconds. If this parameter is omitted or an
%                  empty matrix is passed, then the value of the function
%                  getEOP will be used.
%        clockLoc  An optional 3X1 vector specifying the location of the
%                  clock in WGS-84 ECEF Cartesian [x;y;z] coordinates with
%                  units of meters. Due to relativistic effects, clocks
%                  that are synchronized with respect to TT are not
%                  synchronized with respect to TDB. If this parameter is
%                  omitted or an empty matrix is passed, then a clock at
%                  the center of the Earth is used.
%
%OUTPUTS:   Jul1, Jul2  Matrices of the time as pseudo-Julian dates in TDB
%                       where each row/column corresponds to the values in
%                       the same row/column of the input matrices.
%
%This function just calls the function Cal2TT and then TT2TDB. If deltaT is
%not passed, then getEOP is used to find it.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<8||isempty(clockLoc))
        clockLoc=[0;0;0];
    end

    [Jul1,Jul2]=Cal2TT(year,month,day,hour,minute,second);
    if(nargin<7||isempty(deltaT))
        [JulUTC1,JulUTC2]=TT2UTC(Jul1,Jul2);
        [~,~,~,deltaT]=getEOP(JulUTC1,JulUTC2);
    end

    [Jul1,Jul2]=TT2TDB(Jul1,Jul2,deltaT,clockLoc);
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
