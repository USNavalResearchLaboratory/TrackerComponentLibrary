function [vec,rotMat]=TEME2ITRS(x,Jul1,Jul2,deltaTTUT1,xpyp,LOD)
%%TEME2ITRS Convert from the True Equator Mean Equinox (TEME) of date 
%          coordinate system to the International Terrestrial Reference
%          System (ITRS). The TEME system is non-standard and is generally
%          only used in the Specialized General Perturbations 4 (SGP4)
%          orbit propagation algorithm. Note that the velocity conversion
%          does not include the (small) centrifugal effect of polar
%          motion.
%
%INPUTS: x The NXnumVec collection of vectors in TEME coordinates to
%          convert. N can be 3, or 6. If the vectors are 3D, then they are
%          position. 6D vectors are assumed to be position and velocity,
%          whereby the angular velocity of the Earth's rotation is taken
%          into account using a non-relativistic formula.
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
%     xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
%          including the effects of tides and librations. If this
%          parameter is omitted or an empty matrix is passed the value
%          from the function getEOP will be used.
%      LOD The difference between the length of the day using terrestrial
%          time, international atomic time, or UTC without leap seconds
%          and the length of the day in UT1. This is an instantaneous
%          parameter (in seconds) proportional to the rotation rate of the
%          Earth. This is only needed if more than just position
%          components are being converted.
%
%OUTPUTS: vec A 3XN or 6XN matrix of vectors converted from GCRS
%             coordinates to ITRS coordinates.
%      rotMat The 3X3 rotation matrix used for the conversion of the
%             positions.
%
%The conversion from the TEME to the pseudo-Earth-Fixed (PEF) coordinate
%system is described in [1] and the relationship between the ITRS and the
%PEF is described in [2].
%
%The velocity transformation deals with the instantaneous rotational
%velocity of the Earth using a simple Newtonian velocity addition.
%Basically, the axis of rotation in the Pseudo-Ears-Fixed (PEF) frame is
%the z-axis (The PEF is akin to a less-accurate version of the TIRS).
%The rotation rate in that system is Constants.IERSMeanEarthRotationRate
%adjusted using the Length-of-Day (LOD) Earth Orientation Parameter (EOP).
%Thus, in the PEF, the angular velocity vector is [0;0;omega], where omega
%is the angular velocity accounting for the LOD EOP. Consequently, one
%account for rotation by transforming from the TEME to the PEF,
%subtracting the cross product of Omega with the position in the PEF, and
%then converting to the ITRS. This is a simple Newtonian conversion.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[vec,rotMat]=TEME2ITRS(x,Jul1,Jul2);
%or if more parameters are known,
%[vec,rotMat]=TEME2ITRS(x,Jul1,Jul2,deltaTTUT1,xpyp,LOD);
%
%REFERENCES:
%[1] D. A. Vallado, P. Crawford, R. Hujsak, and T. Kelso, "Implementing
%    the revised SGP4 in STK," in Proceedings of the AGI User Exchange,
%    Washington, DC, 17-18 Oct. 2006, slides. [Online].
%    Available: http://www.agi.com/downloads/events/2006-agi-user-exchange/8_revised_sgp4_vallado2.pdf
%[2] D. A. Vallado, J. H. Seago, and P. K. Seidelmann, "Implementation
%    issues surrounding the new IAU reference systems for astrodynamics,"
%    in Proceedings of the 16th AAS/AIAA Space Flight Mechanics
%    Conference, Tampa, FL, 22-26 Jan. 2006. [Online].
%    Available: http://www.centerforspace.com/downloads/files/pubs/AAS-06-134.pdf
%
%November 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
