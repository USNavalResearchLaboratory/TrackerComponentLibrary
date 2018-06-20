function [xState,errorState]=propagateOrbitSGP4(SGP4Elements,deltaT,TTEpoch1,TTEpoch2,opsMode,gravityModel)
%%PROPAGATEORBITSGP4 Use the SGP4/SDP4 propagator to obtain the position
%                    and velocity of a non-maneuvering Earth-orbiting
%                    satellite after a given time offset. The SGP4
%                    elements are given in the obsolete True Equator Mean
%                    Equinox (TEME) of date coordinate system as is the
%                    result. The propagator is of interest,
%                    because two-line element (TLE) sets, which are sets
%                    of satellite ephemerides published by the U. S. Air
%                    Force are given in terms of SGP4 orbital elements.
%                    The function TLE2SGP4OrbEls  can be used to
%                    get orbital elements for this function from a TLE.
%                    The propagator is actually a combination of the
%                    Simplified General Perturbations 4 (SGP4) dynamic
%                    model for satellites that remain close to the Earth
%                    and the Simplified Deep Space Perturbations 4 (SDP4)
%                    model for satellites that go far enough from the
%                    Earth that the gravitational effects of the Sun and
%                    Moon begin to matter. The propagator only works with
%                    satellites in orbit, not those on an escape
%                    trajectory and is inaccurate with decaying
%                    trajectories. The accuracy of the deep space
%                    propagator, which is used when
%                    2*pi/SGP4Elements(6)>=225, is not very good.
% 
%INPUTS: SGP4Elements A 7X1 or 1X7 vector of SGP4 orbital elements for the
%                     SGP4 propagator that are close to but not the same
%                     as Keplerian orbital elements having the same names.
%                     These are:
%                     1) eccentricity (0-1)
%                     2) inclination  (radians)
%                     3) argument of perigee/periapsis (radians)
%                     4) right ascension of the ascending node (radians)
%                     5) mean anomaly (radians)
%                     6) mean motion (radians per second [TT])
%                     7) BSTAR drag term. This pseudo ballistic
%                       coefficient has units of inverse Earth radii.
%                       Normally, a ballistic coefficient BC is mass per
%                       unit area divided by a drag coefficient. The BSTAR
%                       drag term is a scaled version of the normal
%                       ballistic coefficient. That is,
%                       BSTAR=BC*rho0/2 where rho0=2.461e-5 kg/m^2. The
%                       ballistic coefficient BC is Cd*Area/mass where Cd
%                       is the drag coefficient of the object.
%              deltaT An NX1 or 1XN vector of the time offsets in seconds
%                     (TT) from the epoch time at which one wishes to
%                     determine the state of the target. Note that there
%                     are exactly 86400 TT seconds in a TT Julian day.
%   TTEpoch1, TTEpoch2 The epoch time of the SGP4Elements given in
%                     terrestrial time (TT) as a two-part Julian date (as
%                     a fractional day count). These are only required if
%                     2*pi/SGP4Elements(6)>=225, which is when the deep
%                     space propagator is used. Otherwise, these inputs
%                     can be omitted or empty matrices can be passed. The
%                     added precision of using a split date does not
%                     matter in this instance, since Vallado's SGP4
%                     library routine just uses the summed value.  
%             opsMode An optional parameter specifying the orbital
%                     propagation mode. Possible values are
%                     0 (The default if omitted) Use a model that is
%                       supposed to be similar to what the Air Force Space
%                       Command (AFSPC) uses when they publish orbital
%                       elements.
%                     1 Use a model that is supposed to be a slight
%                       improvement over what the AFSPC uses.
%        gravityModel An optional parameter specifying which gravitational
%                     model should be used. Since the gravitational models
%                     stem from ellipsoidal Earth models, the choice of
%                     the model also affects the radius of the Earth term
%                     used in the algorithm. Possible values are
%                     0 (The default if omitted). Use a gravitational
%                       model based on the WGS-72 reference ellipsoid.
%                       This appears to be what the AFSPC uses in their
%                       TLE sets.
%                     1 Use a gravitational model based on the WGS-84
%                       reference ellipsoid.
%        
%OUTPUTS:      xState The 6XN target state at deltaT offset from the epoch
%                     time. The state consists of position and velocity
%                     components in the order [x;y;z;xDot;yDot;zDot] and
%                     is in the obsolete TEME coordinate system having
%                     units of meters (for position) and meters per
%                     second (for velocity).
%          errorState This indicates whether any errors occurred during
%                     the propagation. Possible values are
%                      0: There is no problem.
%                      1: Mean element problem. Eccentricity >= 1.0 or
%                              eccentricity < -0.001 or semimajor axis
%                              implied by the elements < 0.95.
%                      2: Mean motion less than 0.
%                      3: Partially osculating element problem (these are
%                         derived from the SGP4Elements), eccentricity <
%                         0.0  or  eccentricity > 1.0
%                      4: Semi-latus rectum of the orbit implied by the
%                             elements is < 0.
%                      5: Epoch elements are sub-orbital. This error is
%                         no longer used.
%                      6: Satellite has decayed.
%
%The SGP4 propagator is described in [1] and [2]. The implementation here
%uses the code that Vallado released into the public domain with his
%paper, and which was downloaded from
%http://www.centerforspace.com/downloads/
%as an external library for performing the key propagation step.
%
%Note that this is NOT the official SGP4 orbital propagator used by the
%U.S. Air Force and cannot be assumed to be as reliable or produce
%identical results to the official propagator. Information on obtaining
%the U.S. Air Force's official propagator is given at
%http://www.afspc.af.mil/units/ASDA/
%
%The algorithm can be compiled for use in Matlab  using the
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[xState,errorState]=propagateOrbitSGP4(SGP4Elements,deltaT,TTEpoch1,TTEpoch2,opsMode,gravityModel);
%or as
%[xState,errorState]=propagateOrbitSGP4(SGP4Elements,deltaT,TTEpoch1,TTEpoch2);
%or if no deep-space model is needed, as
%[xState,errorState]=propagateOrbitSGP4(SGP4Elements,deltaT);
%
%REFERENCES:
%[1] F. R. Hoots and R. L. Roehrich, "Spacetrack report no. 3: Models for
%    propagation of NORAD element sets," Department of Defense, Tech.
%    Rep., 31 Dec. 1988. [Online].
%    Available: http://www.amsat.org/amsat/ftp/docs/spacetrk.pdf
%[2] D. A. Vallado, P. Crawford, R. Hujsak, and T. S. Kelso, "Revisiting
%    spacetrack report # 3: Rev 2," in Proceedings of the AIAA/AAS
%    Astrodynamics Specialist Conference and Exhibit, Keystone, CO, 21-24
%    Aug. 2006. [Online].
%    Available: http://celestrak.com/publications/AIAA/2006-6753/AIAA-2006-6753.pdf
%
%December 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
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
