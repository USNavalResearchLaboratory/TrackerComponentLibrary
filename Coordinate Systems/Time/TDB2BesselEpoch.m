function bessEpoch=TDB2BesselEpoch(Jul1,Jul2)
%%TDB2BESSELEPOCH Convert a two-part Julian date in barycentric dynamical
%                  time (TDB) into a Besselian epoch. Besselian years are
%                  related to the apparent position of the Sun and are
%                  generaly considered outdated.
%
%INPUTS: Jul1, Jul2 Matrices of two parts of a Julian date given in TDB.
%                   The units of the date are days. The full date is the
%                   sum of both terms. The date is broken into two parts to
%                   provide more bits of precision. It does not matter how
%                   the date is split. Corresponding elements in each
%                   matrix are times that are converted.
%
%OUTPUTS: bessEpoch The time as a Besselian epoch with the same
%                   dimensionality as the input sets of dates.
%
%The relationship between Julian dates in TDB and Besselian epochs is given
%in [1]. A Besselian epoch is a factional year number denominated in terms
%of a tropical year. For example, 1900.5 is a time near the middle of the
%year 1900. A Besselian year is not the same duration as a Julian year in
%TT.
%
%This is a mex wrapper for the function iauEpb in the International
%Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
%Though the IAU's library does not explicitly say that the Julian date
%must be in TDB, the Besselian epoch does not make sense unless the date is
%provided in TDB as one can ascertain from [1]. 
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%bessEpoch=TDB2BesselEpoch(Jul1,Jul2);
%
%REFERENCES:
%[1] J. H. Lieske, "Precession matrix based on IAU (1976) system of
%    astronomical constants," Astronomy and Astrophysics, vol. 73, no. 3,
%    pp. 282-284, Mar. 1979.
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
