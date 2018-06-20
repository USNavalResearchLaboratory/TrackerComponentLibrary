function [Jul1,Jul2]=TAI2TT(Jul1,Jul2)
%%TAI2TT Convert from international atomic time (TAI) given as a two-part
%        Julian date to terrestrial time (TT), represented as a two-part
%        Julian date.
%
%INPUTS: Jul1, Jul2 Matrices of two parts of a Julian date given in TAI.
%                   The units of the date are days. The full date is the
%                   sum of both terms. The date is broken into two parts
%                   to provide more bits of precision. It does not matter
%                   how the date is split. Corresponding elements in each
%                   matrix are times that are converted.
%
%OUTPUTS: Jul1, Jul2 The time as a Julian date in TT with the same
%                    dimensionalities as the input sets of dates.
%
%This is a mex wrapper for the function iauTaitt in the International
%Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[Jul1,Jul2]=TAI2TT(Jul1,Jul2);
%
%Many temporal coordinate systems standards are compared in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report,
%    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
%    173 pages.
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
