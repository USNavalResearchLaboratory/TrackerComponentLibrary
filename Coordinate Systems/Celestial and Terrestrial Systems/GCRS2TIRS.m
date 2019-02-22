function [vec,rotMat]=GCRS2TIRS(x,Jul1,Jul2,deltaTTUT1,dXdY,LOD)
%%GCRS2TIRS Convert vectors of position and possibly velocity from the
%           Geocentric Celestrial Reference System (GCRS), a type of
%           Earth-Centered Inertial (ECI) coordinate system, to the
%           Terrestrial Intermediate Reference System (TIRS), a type of
%           approximately Earth-Centered Earth-Fixed (ECEF) coordinate
%           system that does not move with polar motion. Adjusting for
%           polar motion would complete the conversion to the
%           International Terrestrial Reference System (ITRS). The
%           conversion of velocity omits the centrifugal effects of the
%           conversion from the GCRS into the Celestial Intermediate
%           Reference System (CIRS). The period of motion of the
%           Celestial Intermediate Pole (CIP) in the GCRS is on the order
%           of 14 months and thus the effect is small.
%
%INPUTS: x The NXnumVec collection of vectors to convert. N can be 3, or
%          6. If the vectors are 3D, then they are position. 6D vectors
%          are assumed to be position and velocity, whereby the angular
%          velocity of the Earth's rotation is taken into account using a
%          non-relativistic formula.
% Jul1, Jul2 Two parts of a Julian date given in terrestrial time (TT).
%          The units of the date are days. The full date is the sum of
%          both terms. The date is broken into two parts to provide more
%          bits of precision. It does not matter how the date is split.
% deltaTTUT1 An optional parameter specifying the difference between TT
%          and UT1 in seconds. This information can be obtained from
%          http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%          or 
%          http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%          If this parameter is omitted or if an empty matrix is passed,
%          then the value provided by the function getEOP will be used
%          instead.
%     dXdY dXdY=[dX;dY] are the celestial pole offsets with respect to
%          the IAU 2006/2000A precession/nutation model in radians If this
%          parameter is omitted, the value from the function getEOP will
%          be used.
%      LOD The difference between the length of the day using terrestrial
%          time, international atomic time, or UTC without leap seconds
%          and the length of the day in UT1. This is an instantaneous
%          parameter (in seconds) proportional to the rotation rate of the
%          Earth. This is only needed if more than just position
%          components are being converted.
%
%OUTPUTS: vec A 3XN or 6XN matrix of vectors converted from GCRS
%             coordinates to TIRS coordinates.
%      rotMat The 3X3 rotation matrix used for the conversion of the
%             positions.
%
%The conversion functions from the International Astronomical Union's
%(IAU) Standard's of Fundamental Astronomy library are put together to get
%the necessary rotation matrix for the position.
%
%The velocity transformation deals with the instantaneous rotational
%velocity of the Earth using a simple Newtonian velocity addition.
%Basically, the axis of rotation in the Terrestrial Intermediate Reference
%System TIRS is the z-axis. The rotation rate in that system is
%Constants.IERSMeanEarthRotationRate adjusted using the Length-of-Day
%(LOD) Earth Orientation Parameter (EOP). Thus, in the TIRS, the angular
%velocity vector is [0;0;omega], where omega is the angular velocity
%accounting for the LOD EOP. Consequently, one accounts for rotation by
%transforming from the GCRS to the TIRS, and subtracting the cross product
%of Omega with the position in the TIRS.
%This is a simple Newtonian conversion.
%
%The algorithm can be compiled for use in Matlab using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[vec,rotMat]=GCRS2TIRS(x,Jul1,Jul2);
%or if more parameters are known,
%[vec,rotMat]=GCRS2TIRS(x,Jul1,Jul2,deltaTTUT1,dXdY,LOD);
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

