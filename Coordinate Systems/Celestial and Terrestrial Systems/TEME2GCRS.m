function [vec,rotMat]=TEME2GCRS(x,Jul1,Jul2,deltaTTUT1,xpyp,dXdY,LOD)
%%TEME2GCRS Convert from the True Equator Mean Equinox (TEME) of date 
%           coordinate system to the Geocentric Celestrial Reference System
%           (GCRS), a type of Earth-Centered Inertial (ECI) coordinate
%           system. The TEME system is non-standard and is generally
%           only used in the Specialized General Perturbations 4 (SGP4)
%           orbit propagation algorithm.
%
%INPUTS: x The NXnumVec collection of vectors in TEME coordinates to
%          convert. N can be 3, or 6. If the vectors are 3D, then they are
%          position. 6D vectors are assumed to be position and velocity,
%          whereby the angular velocity of the Earth's rotation is taken
%          into account using a non-relativistic formula.
% Jul1, Jul2 Two parts of a Julian date given in terrestrial time (TT).
%          The units of the date are days. The full date is the sum of both
%          terms. The date is broken into two parts to provide more bits of
%          precision. It does not matter how the date is split.
% deltaTTUT1 An optional parameter specifying the difference between TT and
%          UT1 in seconds. This information can be obtained from
%          http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%          or 
%          http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%          If this parameter is omitted or if an empty matrix is passed,
%          then the value provided by the function getEOP will be used
%          instead.
%     xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
%          including the effects of tides and librations. If this parameter
%          is omitted or if an empty matrix is passed, the value from the
%          function getEOP will be used.
%     dXdY dXdY=[dX;dY] are the celestial pole offsets with respect to the
%          IAU 2006/2000A precession/nutation model in radians If this
%          parameter is omitted or if an empty matrix is passed, the value
%          from the function getEOP will be used.
%      LOD The difference between the length of the day using terrestrial
%          time, international atomic time, or UTC without leap seconds and
%          the length of the day in UT1. This is an instantaneous parameter
%          (in seconds) proportional to the rotation rate of the Earth.
%          This is only needed if more than just position components are
%          being converted.
%
%OUTPUTS: vec A 3XN or 6XN matrix of vectors converted from TEME
%             coordinates to GCRS coordinates.
%      rotMat The 3X3 rotation matrix used for the conversion of the
%             positions. That is, vec(1:3)=rotMat*x(1:3)
%
%This function just calls the functions TEME2ITRS and ITRS2GCRS, reusing
%the same Earth orientation parameters across calls.
%
%Different celestial coordinate systems are compared in [1].
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report, Naval
%    Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016, 173
%    pages.
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If more EOPs are needed.
if(nargin<7&&size(x,1)==6||nargin<6||isempty(dXdY)||isempty(xpyp)||isempty(deltaTTUT1))
    [JulUTC1,JulUTC2]=TT2UTC(Jul1,Jul2);
    [xpypget,dXdYget,~,deltaTTUT1get,LODget]=getEOP(JulUTC1,JulUTC2);
end

if(nargin<7||isempty(LOD))
   LOD=LODget;
end

if(nargin<6||isempty(dXdY))
    dXdY=dXdYget;
end

if(nargin<5||isempty(xpyp))
    xpyp=xpypget;
end

if(nargin<4||isempty(deltaTTUT1))
    deltaTTUT1=deltaTTUT1get;
end

[x,rotMat1]=TEME2ITRS(x,Jul1,Jul2,deltaTTUT1,xpyp,LOD);
[vec,rotMat2]=ITRS2GCRS(x,Jul1,Jul2,deltaTTUT1,xpyp,dXdY,LOD);

%The combined rotation matrix.
rotMat=rotMat2*rotMat1;
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
