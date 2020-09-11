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
%OUTPUTS: xState The 6XN target state at deltaT offset from the epoch time.
%                The state consists of position and velocity components in
%                the order [x;y;z;xDot;yDot;zDot] and is in the obsolete
%                TEME coordinate system having units of meters (for
%                position) and meters per second (for velocity).
%     errorState This indicates whether any errors occurred during the
%                propagation. Possible values are
%                0: There is no problem.
%                1: Mean element problem. Eccentricity >= 1.0 or
%                   eccentricity < -0.001 or semimajor axis implied by the
%                   elements < 0.95.
%                2: Mean motion less than 0.
%                3: Partially osculating element problem (these are derived
%                   from the SGP4Elements), eccentricity < 0.0  or 
%                   eccentricity > 1.0
%                4: Semi-latus rectum of the orbit implied by the elements
%                   is < 0.
%                5: Epoch elements are sub-orbital. This error is no longer
%                   used.
%                6: Satellite has decayed.
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
%EXAMPLE:
%Given four sets of two-line element sets, determine and plot the
%geometries over the period of one day in ECI coordinates. The satellite
%locations at the start time are marked.
% %GPS BIIF-8  (PRN 03)    
% TLELine1{1}='1 40294U 14068A   15020.29385414  .00000070  00000-0  00000+0 0   540';
% TLELine2{1}='2 40294  54.9733 195.4533 0009937 200.3481 159.6074  2.00552103  1649';
% %MOLNIYA 1-52            
% TLELine1{2}='1 13012U 81123A   15018.94905608 -.00000681  00000-0 -14916+0 0  9748';
% TLELine2{2}='2 13012  64.5459 237.9786 6582041 266.1991 281.9349  2.00764916242351';
% %IRIDIUM 21 [+]          
% TLELine1{3}='1 25778U 99032B   15020.52891112  .00000090  00000-0  24907-4 0   692';
% TLELine2{3}='2 25778  86.3940 231.2123 0002236 103.3219 284.0147 14.34221960822175';
% %TDRS 5 
% TLELine1{4}='1 21639U 91054B   15018.78954554  .00000089  00000-0  00000+0 0  7277';
% TLELine2{4}='2 21639  13.5697  31.2110 0016573 318.7792 244.3708  1.00275730 85939';
% 
% numSat=4;
% SGP4Elements=cell(numSat,1);
% SGP4TTEpoch1=zeros(numSat,1);
% SGP4TTEpoch2=zeros(numSat,1);
% for curTLE=1:numSat
%     [SGP4Elements{curTLE},SGP4TTEpoch1(curTLE),SGP4TTEpoch2(curTLE)]=TLE2SGP4OrbEls(TLELine1{curTLE},TLELine2{curTLE});
% end
% 
% minEpochFrac=min(SGP4TTEpoch2);
% %The constant is the number of seconds in a TT Julian day. We are
% %converting from TT Julian days to seconds.
% SGP4DeltaEpoch=(minEpochFrac-SGP4TTEpoch2)*86400;
% 
% a=Constants.WGS72SemiMajorAxis;%The equatorial radius of the Earth = semi-major axis.
% f=Constants.WGS72Flattening;%The flattening factor fo the Earth.
% b=a*(1-f);%The semi-minor axis.
% 
% %Determine the geometries over a day.
% numPoints=500;
% deltaTOffset=linspace(0,86400,numPoints);%86400 seconds per day.
% 
% xState=zeros(6,numPoints,numSat);
% for curSat=1:numSat
%     for curOffset=1:numPoints
%         deltaT=SGP4DeltaEpoch(curSat)+deltaTOffset(curOffset);
%         %Get the target state in GCRS coordinates at the given time.
%         xTemp=propagateOrbitSGP4(SGP4Elements{curSat},deltaT,SGP4TTEpoch1(curSat),SGP4TTEpoch2(curSat));
%         xState(:,curOffset,curSat)=TEME2GCRS(xTemp,SGP4TTEpoch1(curSat),SGP4TTEpoch2(curSat)+deltaT/86400);
%     end
% end
% 
% %Plot the geometries over the period of one day.
% figure(1)
% clf
% hold on
% [x,y,z]=ellipsoid(0,0,0,a,a,b,20);
% surf(x,y,z)
% colormap(gray);
% 
% %Orbits
% plot3(xState(1,:,1),xState(2,:,1),xState(3,:,1),'-b','linewidth',4)
% plot3(xState(1,:,2),xState(2,:,2),xState(3,:,2),'-r','linewidth',4)
% plot3(xState(1,:,3),xState(2,:,3),xState(3,:,3),'-g','linewidth',4)
% plot3(xState(1,:,4),xState(2,:,4),xState(3,:,4),'-k','linewidth',4)
% 
% %Satellite Initial Locations
% scatter3(xState(1,1,1),xState(2,1,1),xState(3,1,1),100,'sb','filled','MarkerEdgeColor','b','linewidth',3)
% scatter3(xState(1,1,2),xState(2,1,2),xState(3,1,2),100,'or','filled','MarkerEdgeColor','r','linewidth',3)
% scatter3(xState(1,1,3),xState(2,1,3),xState(3,1,3),100,'og','filled','MarkerEdgeColor','g','linewidth',3)
% scatter3(xState(1,1,4),xState(2,1,4),xState(3,1,4),100,'ok','filled','MarkerEdgeColor','k','linewidth',3)
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
