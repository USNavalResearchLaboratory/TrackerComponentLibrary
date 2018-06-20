function [xHelio,xBary,exitFlag]=approxSolarSysVec(Jul1,Jul2,solarBody)
%%APPROXSOLARSYSVEC Get the approximate position and velocity vector of
%                 the Earth or of the planets in coordinates aligned with
%                 the International Celestial Reference System (ICRS).
%                 The position for everything is heliocentric, and for the
%                 Earth an additional barycentric (Like the Barycentric
%                 Celestial Reference System [BCRS], except TDB is used
%                 and not TCB) vector is returned. Use the function
%                 solarBodyVec if higher-precision with respect to
%                 observers near the Earth is desired.
%
%INPUTS: Jul1,Jul2 Two parts of a Julian date given in Barycentric
%                  dynamical time (TDB). However, terrestrial time (TT)
%                  can generally be substituted as the two timescales are
%                  relatively close and the approximation is not very high
%                  fidelity. The units of the date are days. NX1 or 1XN 
%                  vectors can be passed if one wishes to find the 
%                  positions and velocities at multiple times. The full 
%                  date is the sum of both terms. The date is broken into
%                  two parts to provide more bits of precision. 
%                  It does not matter how the date is split. The 
%                  valid range of dates is from 1900AD-2100AD for the 
%                  Earth (solarBody=0) and 1000AD- 3000AD for everything
%                  else.
%        solarBody A parameter specifying the solar body. Possible values
%                  are:
%                  0 (The default if omitted) The Earth.
%                  1 Mercury
%                  2 Venus
%                  3 The Earth-Moon barycenter
%                  4 Mars
%                  5 Jupiter
%                  6 Saturn
%                  7 Uranus
%                  8 Neptune
%
%OUTPUTS: xHelio  A 6XN state vector of the body where the coordinate
%                 axes are aligned with the ICRS and the origin is the
%                 Sun. xHelio(1:3,i) is the position for the ith date and
%                 xHelio(4:6,i) is the velocity. If an error occurred for
%                 a particular date, zeros are returned for that date.
%                 The units are meters and meters per second TDB (which is
%                 close to TT).
%         xBary   A 6XN state vector of the body where the coordinate
%                 axes are aligned with the ICRS and the origin is the
%                 solar system barycenter. This is essentially the BCRS.
%                 xBary(1:3,i) is the position for the ith date and
%                 xBary(4:6,i) is the velocity. If an error occurred for
%                 a particular date, zeros are returned for that date.
%                 The units are meters and meters per second TDB (which is
%                 close to TT). xBary is only returned if solarBody=0.
%                 Otherwise, an empty matrix is returned.
%        exitFlag A NX1 set of numbers indicating the exit status.
%                 Possible values for each date are
%                 0 = OK
%                 1 = Date outside of valid range.
%                 2 = Algorithm failed to converge.
%
%This function is just a wrapper for the iauEpv00 function in the
%International Astronomical Union's (IAU) Standard's of Fundamental
%Astronomy library.
%
%The algorithm can be compiled for use in Matlab  using the 
%CompileCLibraries function.
%
%The algorithm is run in Matlab using the command format
%[xHelio,xBary,exitFlag]=approxSolarSysVec(Jul1,Jul2,solarBody);
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
