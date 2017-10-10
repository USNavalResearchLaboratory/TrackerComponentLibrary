function [Jul1,Jul2]=TAI2TCB(Jul1,Jul2,deltaT,clockLoc)
%%TAI2TCB Convert from international atomic time (TAI) to barycentric
%         coordinate time (TCB) to an accuracy of nanoseconds (if deltaT
%         is accurate) using the routines from the International
%         Astronomical Union's library that do not require external
%         ephemeris data.
%
%INPUTS: Jul1,Jul2 Two parts of a Julian date given in TT. The units of
%                  the date are days. The full date is the sum of both
%                  terms. The date is broken into two parts to provide
%                  more bits of precision. It does not matter how the date
%                  is split.
%           deltaT An optional parameter specifying the offset between TT
%                  and UT1 in seconds. If this parameter is omitted or an
%                  empty matrix is passed, then the value of the function
%                  deltaTTUT1 will be used.
%         clockLoc An optional 3X1 vector specifying the location of the
%                  clock in WGS-84 ECEF Cartesian [x;y;z] coordinates with
%                  units of meters. Due to relativistic effects, clocks
%                  that are synchronized with respect to TAI are not
%                  synchronized with respect to TCB. If this parameter is
%                  omitted, then a clock at the center of the Earth is
%                  used and the precision declines to microseconds.
%        
%OUTPUTS:Jul1,Jul2 Two parts of a Julian date given in TCB.
%
%This function converts the date to a terrestrial time and then calls
%TT2TCB.
%
%Many temporal coordinate systems standards are compared in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report, Naval
%    Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016, 173
%    pages.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[Jul1,Jul2]=TAI2TT(Jul1,Jul2);

switch(nargin)
    case 2
        [Jul1,Jul2]=TT2TCB(Jul1,Jul2);
    case 3
        [Jul1,Jul2]=TT2TCB(Jul1,Jul2,deltaT);
    otherwise
        [Jul1,Jul2]=TT2TCB(Jul1,Jul2,deltaT,clockLoc);
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
