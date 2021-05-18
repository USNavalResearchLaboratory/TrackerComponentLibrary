function [Jul1,Jul2]=JulEpoch2JulDate(epochJ)
%%JULEPOCH2JULDATE Convert a Julian epoch to a Julian Date in the same
%                  timescale. For example, an epoch given in terrestrial
%                  time (TT) will be a Julian Date in TT. The timescale
%                  should be uniform (i.e. not UTC).
%
%INPUTS: epochJ A matrix of Julian epochs as a fractional year in a
%               uniform timescale.
%
%OUTPUTS: Jul1, Jul2 Two parts of a Julian date given in the same
%                    timescale as the epoch The units of the date are
%                    days. The full date is the sum of both terms. The
%                    entries in the matrices correspond to the entries in
%                    epochJ.
%
%A Julian epoch is a factional year number denominated in terms of a
%year of exactly 365.25 days. For example 2000.5 is a Julian epoch
%in the middle of the year 2000 in a particular timescale.
%
%This is a wrapper for the function iauEpj2jd in the International
%Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
%
%The algorithm can be compiled for use in Matlab  using the
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[Jul1,Jul2]=JulEpoch2JulDate(epochJ);
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
