function catData=changeEpoch(catData,Jul1,Jul2,Jul3,Jul4)
%%CHANGEEPOCH Given star catalog data for a particular epoch, change it to
%             a new epoch. This function is useful for transforming data
%             in star catalogs to the J2000.0 epoch so that the function
%             starCat2Obs can be used to determine what an observer on the
%             Earth should see.
%
%INPUTS: catData Cat data is a matrix of stars (one per row) that are to
%           be converted, where each row has the following format:
%  catData(:,1) RArad Right Ascension in (rad) ICRS at the J2000.0 epoch
%  catData(:,2) DErad Declination in (rad) ICRS at the J2000.0 epoch
%  catData(:,3) Plx   Parallax (rad)
%  catData(:,4) pmRA  Proper motion in Right Ascension (rad/yr) in the form
%                     dRA/dt and not cos(Dec)*dRA/dt.
%  catData(:,5) pmDE  Proper motion in Declination (rad/yr)
%  catData(:,6) vRad  Radial velocity of the star in meters per second with
%                     a positive sign meaning that a star is receding
%                    (going away from) Earth.
% Jul1,Jul2 Two parts of a Julian date given in TDB specifying the date of
%           the original epoch. The units of the date are days. The full
%           date is the sum of both terms. The date is broken into two
%           parts to provide more bits of precision. It does not matter how
%           the date is split. If these parameters are omitted or if empty
%           matrices are passed, then a source epoch of Jul1=2448349;
%           Jul2=0.062500019022333; is used. These values are what one
%           obtains by putting the terrestrial time 2448349.0625, the epoch
%           of the Hipparcos catalog, though the function TT2TDB.
% Jul3,Jul4 Two parts of a Julian date for the epoch into which the
%           parameters should be converted in TDB. If these parameters are
%           omitted, then a terrestrial time of Jul1=2451545.0; Jul2=0;
%           will be used, corresponding to the J2000.0 epoch.
%
%OUTPUTS: catData The same as the input catData, except the epoch has
%                 been changed.
%
%This function is pretty much a wrapper for the iauPmsafe that is part of
%the International Astronomical Union's Standards of Fundamental Astronomy
%library.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%catData=changeEpoch(catData,Jul1,Jul2,Jul3,Jul4);
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
