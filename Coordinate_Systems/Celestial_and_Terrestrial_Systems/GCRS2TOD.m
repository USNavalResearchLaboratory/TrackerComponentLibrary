function [xRot,rotMat]=GCRS2TOD(xVec,TT1,TT2,dXdY)
%%GCRS2TOD Rotate a vector from the geocentric celestial reference system
%          (GCRS) to the true equator and equinox of date coordinate
%          system (TOD) using the IAU 2006/2000A model, where a location
%          is known as an "apparent place". The transformation is
%          performed by removing the precession, nutation, and frame bias.
%
%INPUTS:  xVec The 3XN matrix of N 3X1 Cartesian vectors that are to be
%              rotated from the GCRS into the TOD.
% TT1, TT2 Jul1,Jul2 Two parts of a Julian date given in TT. The units of
%              the date are days. The full date is the sum of both terms.
%              The date is broken into two parts to provide more bits of
%              precision. It does not matter how the date is split.
%         dXdY dXdY=[dX;dY] are the celestial pole offsets with respect to
%              the IAU 2006/2000A precession/nutation model in radians If
%              this parameter is omitted or an empty matrix is passed, the
%              value from the function getEOP will be used.
%
%OUTPUTS: xRot The 3XN matrix of the N 3X1 input vector rotated into the
%              TOD.
%       rotMat The 3X3 rotation matrix such that
%              xRot(:,i)=rotMat*xVec(:,i).
%
%This uses functions in the International Astronomical Union's (IAU)
%Standard's of Fundamental Astronomy (SOFA) library to obtain the product
%of the nutation and precession rotation matrices  and the frame rotation
%bias matrix. One goes from GCRS to TOD by applying a frame bias and then
%precession and a nutation. The rotations are discussed in the
%documentation for the SOFA library as well as in [1] among other sources.
%
%The correction for using dXdY is the most accurate formula in [2].
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[xRot,rotMat]=GCRS2TOD(xVec,TT1,TT2);
%or
%[xRot,rotMat]=GCRS2TOD(xVec,TT1,TT2,dXdY);
%
%Different celestial coordinate systems are compared in [3].
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%[2] G. H. Kaplan, "Celestial pole offsets: Conversion from (dx,dy) to
%   (dpsi,deps)," U.S. Naval Observatory, Tech. Rep., May 2005. [Online].
%   Available: http://aa.usno.navy.mil/publications/reports/dXdY_to_dpsideps.pdf
%[3] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
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
