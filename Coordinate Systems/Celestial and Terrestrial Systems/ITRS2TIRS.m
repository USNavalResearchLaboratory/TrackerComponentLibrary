function [vITRS,rotMat]=ITRS2TIRS(vTIRS,TT1,TT2,xpyp)
%%ITRS2TIRS Rotate vectors from the International terrestrial Reference
%            System (ITRS) into the  Terrestrial Intermediate Reference
%            System (TIRS). The ITRS is essentially the WGS-84 coordinate
%            system: it defines locations with respect to the crust of a
%            non-rotating Earth, where the z axis passes through a fixed
%            point on the surface. On the other hand, the TIRS is nearly
%            the same except the z axis is the axis of rotation of the
%            Earth, which slowly varies over time. Note that the velocity
%            conversion does not include the (small) centrifugal effect of
%            polar motion.
%
%INPUTS: x The NXnumVec collection of vectors in TIRS coordinates to
%          convert (units do not matter). N can be 3, or 6. If the vectors
%          are 3D, then they are position. 6D vectors are assumed to be
%          position and velocity. Since the TIRS and ITRS co-rotate, there
%          is no Coriolis effect to add. Also, the accelerations due to the
%          wobble of the rotation axis over time are not considered. These
%          accelerations are very small. Thus, the function just rotates
%          both halves of the vector.
% TT1, TT2 Two parts of a Julian date given in terrestrial time (TT).
%          The units of the date are days. The full date is the sum of
%          both terms. The date is broken into two parts to provide more
%          bits of precision. It does not matter how the date is split.
%     xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
%          including the effects of tides and librations. If this
%          parameter is omitted or if an empty matrix is passed, the value
%          from the function getEOP will be used.
%
%OUTPUTS: vITRS The NXnumVec vector of values of x rotated from the ITRS
%               into the TIRS.
%        rotMat The 3X3 rotation matrix used to rotate vectors from the
%               ITRS into the TIRS.
%
%The conversion functions from the International Astronomical Union's
%(IAU) Standard's of Fundamental Astronomy library are put together to get
%the necessary rotation matrix.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[vITRS,rotMat]=ITRS2TIRS(vTIRS,TT1,TT2)
%or if more parameters are known, using the format
%[vITRS,rotMat]=ITRS2TIRS(vTIRS,TT1,TT2,xpyp);
%
%Different celestial coordinate systems are compared in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report,
%    Naval Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016,
%    173 pages.
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
