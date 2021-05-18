function [Jul1,Jul2]=TDB2TT(Jul1,Jul2,deltaTTUT1,clockLoc)
%%TDB2TT Convert from  barycentric dynamical time (TDB) to terrestrial
%        time (TT) to an accuracy of nanoseconds (if deltaT is accurate)
%        using the routines from the International Astronomical Union's
%        library that do not require external ephemeris data.
%
%INPUTS: Jul1,Jul2 Two parts of a Julian date given in TT. The units of
%                  the date are days. The full date is the sum of both
%                  terms. The date is broken into two parts to provide
%                  more bits of precision. It does not matter how the date
%                  is split.
%       deltaTTUT1 An optional parameter specifying the offset between TT
%                  and UT1 in seconds. If this parameter is omitted or an
%                  empty matrix is passed, then the value of the function
%                  getEOP will be used. Since getEOP takes UTC (which
%                  comes from UT1) and we have only TDB (which is close),
%                  a few iterations are performed to try to get the
%                  correct result.
%         clockLoc An optional 3X1 vector specifying the location of the
%                  clock in the Terrestrial Intermediate Reference System
%                  (TIRS), though it would not make much of a difference
%                  if the International Terrestrial Reference System
%                  (ITRS) were used. The units are meters. Due to
%                  relativistic effects, clocks that are synchronized with
%                  respect to TT are not synchronized with respect to TDB.
%                  If this parameter is omitted, then a clock at the
%                  center of the Earth is used.
%        
%OUTPUTS:Jul1,Jul2 Two parts of a Julian date given in TDB.
%
%This function relies on a number of functions in the International
%Astronomical Union's Standards of Fundamental Astronomy library.
% 
%The main implementation issue is that if deltaT is not provided (and
%generally, one probably would not be expected to know it as it is not
%tabulated in TDB), one has to iterate using the getEOP function to get
%the correct offset to use. However, even if the deltaT parameter is
%given, one still has to iterate, because values in UT1 are required and
%the iauTtut1 function requires UT1.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[Jul1,Jul2]=TDB2TT(Jul1,Jul2,deltaTTUT1,clockLoc);
%or
%[Jul1,Jul2]=TDB2TT(Jul1,Jul2);
%
%Many temporal coordinate systems standards are compared in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report,
%    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
%    173 pages.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
