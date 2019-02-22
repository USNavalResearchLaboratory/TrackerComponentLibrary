function bessEpoch=TT2BesselEpoch(Jul1,Jul2,deltaTTUT1,clockLoc)
%%TT2BESSEPOCH Convert from terrestrial time (TT) into a Besselian epoch.
%              Besselian years are related to the apparent position of the
%              Sun and are generaly considered outdated.
%
%INPUTS: Jul1, Jul2 Matrices of two parts of a Julian date given in TT.
%                   The units of the date are days. The full date is the
%                   sum of both terms. The date is broken into two parts to
%                   provide more bits of precision. It does not matter how
%                   the date is split. Corresponding elements in each
%                   matrix are times that are converted.
%        deltaTTUT1 An optional parameter specifying the offset between TT
%                   and UT1 in seconds. If this parameter is omitted or an
%                   empty matrix is passed, then the value of the function
%                   getEOP will be used.
%          clockLoc An optional 3X1 vector specifying the location of the
%                   clock in WGS-84 ECEF Cartesian [x;y;z] coordinates with
%                   units of meters. Due to relativistic effects, clocks
%                   that are synchronized with respect to TT are not
%                   synchronized with respect to TDB. If this parameter is
%                   omitted, then a clock at the center of the Earth is
%                   used.
%
%OUTPUTS: bessEpoch The time as a Besselian epoch with the same
%                   dimensionality as the input sets of dates.
%
%This function just calls TDB2TT followed by TDB2BesselEpoch.
%
%July 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    deltaTTUT1=[];
end

if(nargin<4)
   clockLoc=[]; 
end

[Jul1,Jul2]=TT2TDB(Jul1,Jul2,deltaTTUT1,clockLoc);
bessEpoch=TDB2BesselEpoch(Jul1,Jul2);

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
