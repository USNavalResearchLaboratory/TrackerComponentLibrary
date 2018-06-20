function xG=ICRS2G(x,Jul1,Jul2)
%%ICRS2G Rotate vectors from the orientation of the International
%        Celestial Reference System (ICRS) to that of the International
%        Astronomical Union's (IAU's) 1958 system of galactic coordinates.
%        This function works with Cartesian vectors and with two spherical
%        angles (azimuth and elevation) without range. The ICRS is aligned
%        with the GCRS and BCRS.
%
%INPUTS: x The NXnumVec collection of vectors to convert. N can be 2, or
%          3. If the vectors are 2D, then they are assumed to be azimuth
%          and elevation in radians. 3D vectors are assumed to be
%          Cartesian position.
%
%OUTPUTS: xG The vectors rotated into the galactic coordinate system. If
%            the input was 2D azimuth and elevation, the output will be
%            the same. If the input was Cartesian, then the output will be
%            Cartesian.
%
%This is mostly a wrapper for the function iauIcrs2g in the International
%Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
%
%The original coordinate system is defined in [1]. However, as documented
%in the IAU's function iauIcrs2g, the galactic coordinate system has been
%redefined in terms of more accurate modern data.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%xG=ICRS2G(x,Jul1,Jul2);
%
%REFERENCES:
%[1] A. Blaauw, C. S. Gum, J. L. Pawsey, and G. Westerhout, "The new
%    I.A.U. system of galactic coordinates (1958 revision)," Monthly Notes
%    of the Royal Astronomical Society, vol. 121, no. 2, pp. 123?131,
%    1960.
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
