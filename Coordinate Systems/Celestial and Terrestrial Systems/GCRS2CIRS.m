function [vec,rotMat]=GCRS2CIRS(x,Jul1,Jul2,dXdY)
%%GCRS2CIRS Convert vectors of position and possibly velocity from the
%           Geocentric Celestrial Reference System (GCRS), a type of
%           Earth-Centered Inertial (ECI) coordinate system, to the
%           Celestial Intermediate Reference System (CIRS), which is
%           offset by the rotation of the Celestial Intermediate Pole
%           (CIP). The velocity conversion omits the centrifugal effects
%           of the CIP motion, which have a period on the order of 14
%           months and are thus small.
%
%INPUTS: x The NXnumVec collection of vectors to convert. N can be 3, or
%          6. If the vectors are 3D, then are position. 6D vectors are
%          assumed to be position and velocity.
% Jul1, Jul2 Two parts of a Julian date given in terrestrial time (TT).
%          The units of the date are days. The full date is the sum of
%          both terms. The date is broken into two parts to provide more
%          bits of precision. It does not matter how the date is split.
%     dXdY dXdY=[dX;dY] are the celestial pole offsets with respect to
%          the IAU 2006/2000A precession/nutation model in radians If this
%          parameter is omitted, the value from the function getEOP will
%          be used.
%
%OUTPUTS: vec A 3XN or 6XN matrix of vectors converted from GCRS
%             coordinates to CIRS coordinates.
%      rotMat The 3X3 rotation matrix used for the rotation of the
%             positions and velocities.
%
%The conversion functions from the International Astronomical Union's
%(IAU) Standard's of Fundamental Astronomy library are put together to get
%the necessary rotation matrix for the position.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[vec,rotMat]=GCRS2CIRS(x,Jul1,Jul2);
%or if more parameters are known,
%[vec,rotMat]=GCRS2CIRS(x,Jul1,Jul2,dXdY);
%
%Different celestial coordinate systems are compared in [1].
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
