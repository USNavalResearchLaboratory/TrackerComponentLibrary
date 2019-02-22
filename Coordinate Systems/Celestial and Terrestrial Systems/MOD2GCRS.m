function [xRot,rotMat]=MOD2GCRS(xVec,TT1,TT2)
%%MOD2GCRS Rotate a vector from the mean of date (MOD) coordinate system,
%          which is the coordinate system using the mean equinox and
%          ecliptic of date, to the geocentric celestial reference system
%          (GCRS) using the IAU 2006/2000A model. The transformation is
%          performed by removing the precession and frame bias.
%
%INPUTS: xVec The 3XN matrix of N 3X1 Cartesian vectors that are to be
%             rotated from the MOD coordinate system into the GCRS
%             coordinate system.
%    TT1, TT2 Jul1,Jul2 Two parts of a Julian date given in TT. The units
%             of the date are days. The full date is the sum of both
%             terms. The date is broken into two parts to provide more
%             bits of precision. It does not matter how the date is split.
%
%OUTPUTS: xRot The 3XN matrix of the N 3X1 input vector rotated into the
%               GCRS coordinate system.
%        rotMat The 3X3 rotation matrix such that
%               xRot(:,i)=rotMat*xVec(:,i).
%
%This uses functions in the International Astronomical Union's (IAU)
%Standard's of Fundamental Astronomy (SOFA) library to obtain the product
%of the precession rotation matrix and the frame rotation bias matrix. One
%goes from GCRS to mean of date by applying a frame bias and then
%precession. Thus this function removes those rotations. The rotations are
%discussed in the documentation for the SOFA library as well as in
%G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%Rotation and Reference Systems Service Std. 36, 2010.
%among other sources.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[xRot,rotMat]=MOD2GCRS(xVec,TT1,TT2);
%
%Different celestial coordinate systems are compared in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report,
%    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
%    173 pages.
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
