function eclipVec=ICRS2Eliptic(vec,TT1,TT2,method)
%%ICRS2ECLIPTIC Convert a location vector from the International
%               Celestial Reference System (ICRS) to eliptic coordinates
%               either using the IAU 2006 precession model or the Vondrak
%               400 millennia precession model.
%
%INPUTS: x The NXnumVec collection of vectors in the ICRS to convert. N
%          can be 2, or 3. If the vectors are 2D, then they are assumed to
%          be azimuth and elevation in radians. 3D vectors are assumed to
%          be Cartesian position.
% Jul1, Jul2 Two parts of a Julian date given in terrestrial time (TT).
%          The units of the date are days. The full date is the sum of
%          both terms. The date is broken into two parts to provide more
%          bits of precision. It does not matter how the date is split.
%   method An optional parameter specifying which algorithm is to be used.
%          Possible values are
%          0 (The default if omitted or an empty matrix is passed) Use the
%            IAU 2006 precession model.
%          1 Use the long-term (Vondrak) precession model.
%
%OUTPUTS: xG The vectors rotated into the ecliptic coordinate system. If
%            the input was 2D azimuth and elevation, the output will be
%            the same. If the input was Cartesian, then the output will be
%            Cartesian.
%
%This function is a Matlab interface for the relevant functions in the
%International Astronomical Union's (IAU) Standard's of Fundamental
%Astronomy library.
%
%The ecliptic is defined in the IERS Conventions [1] to be the "the
%plane perpendicular to the mean heliocentric orbital angular momentum
%vector of the Earth-Moon barycentre in the BCRS".
%
%The algorithm can be compiled for use in Matlab  using the
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%eclipVec=ICRS2Eliptic(vec,TT1,TT2,method);
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
