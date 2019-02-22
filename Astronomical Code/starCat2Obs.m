function [zObs,uObs]=starCat2Obs(catData,Jul1,Jul2,zObs,P,T,R,wl,dut1,xpyp)
%%STARCAT2OBS Convert data for the location of stars as typically
%             supplied by a star catalog, such as the Hipparcos catalog,
%             to local ENU observed coordinates at the receiver,
%             including corrections for parallax, aberation,
%             gravitational deflection by the sun and a low-fidelity
%             refraction model. The catalog data is assumed to be at the
%             J2000 epoch.
%
%catData   catData is a matrix of stars (one per row) that are to be
%          converted, where each row has the following format:
%catData(:,1) RArad    Right Ascension in (rad) ICRS at the J2000.0 epoch 
%catData(:,2) DErad    Declination in (rad) ICRS at the J2000.0 epoch
%catData(:,3) Plx      Parallax (rad)
%catData(:,4) pmRA     Proper motion in Right Ascension (rad/yr) in the
%                      form dRA/dt and not cos(Dec)*dRA/dt.
%catData(:,5) pmDE     Proper motion in Declination (rad/yr)
%catData(:,6) vRad     Radial velocity of the star in meters per second
%                      with a positive sign meaning that a star is
%                      receding (going away from) Earth.
% Jul1,Jul2 Two parts of a pseudo-Julian date given in UTC for when the
%          observation is made. The units of the date are days. The full
%          date is the sum of both terms. The date is broken into two
%          parts to provide more bits of precision. It does not matter how
%          the date is split.
%     zObs zObs=[lat;lon;h], the longitude, geodetic latitude and height
%          above the reference ellipsoid of the observer using the WGS-84
%          reference ellipsoid. The units of lat and lon are radians and
%          the height is in meters. East and North are the positive
%          directions.
%        R The relative humidity at the observer (between 0 and 1). If
%          this parameter is omitted or an empty matrix is passed, then
%          Constants.standardRelHumid is used.
%        P The atmospheric pressure at the observer in Pascals (N/m^2). If
%          this parameter is omitted or an empty matrix is passed, then
%          Constants.standardAtmosphericPressure is used.
%        T The air temperature at the observer in degrees Kelvin. If this
%          parameter is omitted or an empty matrix is passed, then
%          Constants.standardTemp is used.
%       wl The wavelength at which the observation is made in units of
%          meters. If this parameter is omitted or an empty matrix is
%          passed, then a wavelength of 0.574 micrometers is used, which
%          is in the visible spectrum (a rather yellow color).
%   deltaT The difference between UTC and UT1 in seconds. This
%          information can be obtained from
%          http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%          or 
%          http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%          If this parameter is omitted or if an empty matrix is passed,
%          then the value provided by the function getEOP will be used
%          instead.
%     xpyp xpyp=[xp,yp], the polar motion coordinates of the respect
%          to the International Terrestrial Reference System. As
%          described in Section 5.1 of the IERS Conventions 2010,
%          values are published by the IERS and should have been
%          updated to account for the additional temporal effects of
%          ocean tides and librations. If this parameter is omitted,
%          the value provided by the function getEOP will be used instead.
%
%OUTPUTS: zSpher For N stars, this is a 2XN matrix with each colum being
%                the location of a star in [azimuth;elevation] in radians
%                taken with respect to a local East-North-Up coordinate
%                system defined on the WGS-84 ellipsoid. Azimuth is
%                measured in radians North of East.
%           uObs For N stars, this is a 3XN matrix of unit vectors in
%                WGS-84 ENU coordinates pointing toward the stars. 
%
%This is a mex wrapper for the function iauAtco13 in the International
%Astronomical Union's (IAU) Standard's of Fundamental Astronomy library.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[zObs,uObs]=starCat2Obs(catData,Jul1,Jul2,zObs,P,T,R,wl,dut1,xpyp);
%or
%[zObs,uObs]=starCat2Obs(catData,Jul1,Jul2,zObs);
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
