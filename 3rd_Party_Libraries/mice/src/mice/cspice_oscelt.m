%-Abstract
%
%   CSPICE_OSCELT calculates the set of osculating conic
%   orbital elements corresponding to the state 6-vector
%   (position, velocity) of a body at an epoch.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY,
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING,
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      state  a double precision 6-vector or 6xN array of
%             states of the body at some epoch. Components
%             are x, y, z, dx/dt, dy/dt, dz/dt. `state' must be
%             expressed relative to an inertial reference frame.
%             Units are km and km/sec.
%
%      et     the double precision scalar or 1XN-vector of ephemeris
%             time epochs corresponding to each 'state' in ephemeris
%             seconds past J2000
%
%      mu     the gravitational parameter of the primary
%             body for 'state'
%
%   the call:
%
%      elts = cspice_oscelt( state, et, mu )
%
%   returns:
%
%      elts   a double precision 8-vector or 8xN array containing
%             the equivalent conic elements describing the orbit
%             of the body around its primary. The elements are,
%             in order:
%
%                 elts(1)  contains rp, perifocal distance.
%                 elts(2)  contains ecc, eccentricity.
%                 elts(3)  contains inc, inclination.
%                 elts(4)  contains lnode, longitude of the ascending node.
%                 elts(5)  contains argp, argument of periapsis.
%                 elts(6)  contains m0, mean anomaly at epoch.
%                 elts(7)  contains t0, epoch.
%                 elts(8)  contains mu, gravitational parameter.
%
%             The epoch of the elements is the epoch of the input
%             state. Units are km, rad, rad/sec. The same elements
%             are used to describe all three types (elliptic,
%             hyperbolic, and parabolic) of conic orbit
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%
%      Example(1):
%
%      %
%      % Determine the osculating elements of the moon wrt the
%      % Earth at some arbitrary time in the J2000 inertial frame.
%      %
%      % Load the meta kernel listing the needed SPK, PCK, LSK
%      % kernels.
%      %
%      cspice_furnsh('standard.tm' )
%
%      %
%      % Convert the time string to ephemeris time
%      %
%      et = cspice_str2et( 'Dec 25, 2007' );
%
%      %
%      % Make the cspice_spkezr call to retrieve the state of the
%      % moon wrt the Earth in J2000.
%      %
%      [state, ltime] = cspice_spkezr( 'Moon', et, 'J2000', 'LT+S', 'EARTH' );
%
%      %
%      % cspice_oscelt requires body mass information, so load a
%      % mass PCK kernel.
%      %
%      cspice_furnsh( '/kernels/gen/pck/masses3.tpc' )
%
%      %
%      % Read the gravitational parameter for Earth.
%      %
%      mu = cspice_bodvrd( 'EARTH', 'GM', 1 );
%
%      %
%      % make the cspice_oscelt call to convert the state 6-vector
%      % to the elts 8-vector. Note: the  cspice_bodvrd returns
%      % data as arrays, so to access the gravitational parameter
%      % (the only value in the array), we use mu(1).
%      %
%      elts = cspice_oscelt( state, et, mu(1) );
%
%      %
%      % Output the elts vector in a column format.
%      %
%      txt = sprintf( '%24.8f\n', elts );
%      disp( txt)
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%         360956.91440370
%              0.07820299
%              0.48717811
%              6.18584105
%              1.28603872
%              0.55386000
%      251812864.18370920
%         398600.44800000
%
%      Example(2):
%
%      %
%      % Calculate the history of the Moon's orbit plane
%      % inclination with respect to the Earth in the
%      % J2000 frame at intervals of one month for a
%      % time interval of 14 years.
%      %
%      % Load the needed SPK , PCK and LSK kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % cspice_oscelt also requires mass information, so load
%      % a mass PCK.
%      %
%      cspice_furnsh( '/kernels/gen/pck/masses3.tpc' )
%      mu = cspice_bodvrd( 'EARTH', 'GM', 1 );
%
%      %
%      % The start epoch.
%      %
%      et0 = cspice_str2et( 'Jan 1, 2000 12:00:00' );
%
%      %
%      % A step of one month - in seconds.
%      %
%      step = 30. * cspice_spd;
%
%      %
%      % Define an array of ephemeris times, covering,
%      % 14 years of months in steps of one month starting
%      % approximately Feb 1, 2000.
%      %
%      et = [0: (14*12) - 1]*step + et0;
%
%      % Retrieve the state; convert to osculating elements.
%      %
%      [state,ltime] = cspice_spkezr( 'Moon', et, 'J2000', 'LT+S', 'EARTH');
%      elts          = cspice_oscelt( state, et, mu(1) );
%
%      elts(3,:) = [ elts(3,:) * cspice_dpr ];
%
%      %
%      % Convert the ephemeris time of the state lookup to
%      % calendar UTC, then print the calendar string and the
%      % inclination in degrees of the Moon wrt Earth at the
%      % time.
%      %
%      utcstr = cspice_et2utc( et, 'C', 3, );
%
%      %
%      % Convert the angular measures to degrees.
%      %
%
%      %
%      % Output the epoch and corresponding inclination.
%      %
%      for n=1:14*12
%         fprintf( '%s %12.6f\n', utcstr(n,:), elts(3,n) );
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      ... a partial output ...
%
%      2012 NOV 23 12:00:00.001    20.903479
%      2012 DEC 23 12:00:00.000    20.902973
%      2013 JAN 22 11:59:59.999    20.802204
%      2013 FEB 21 11:59:59.999    20.565404
%      2013 MAR 23 11:59:59.998    20.309740
%      2013 APR 22 11:59:59.998    20.171117
%      2013 MAY 22 11:59:59.999    20.162453
%      2013 JUN 21 12:00:00.000    20.173366
%      2013 JUL 21 12:00:00.000    20.082464
%      2013 AUG 20 12:00:00.001    19.867300
%      2013 SEP 19 12:00:00.002    19.628911
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine oscelt_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   conic elements from state
%   osculating elements from state
%   convert state to osculating elements
%
%-&

function [elts] = cspice_oscelt( state, et, mu )

   switch nargin
      case 3

         state = zzmice_dp(state);
         et    = zzmice_dp(et);
         mu    = zzmice_dp(mu);

      otherwise

         error ( 'Usage: [_elts(8)_] = cspice_oscelt( _state(6)_, _et_, mu )' )

   end

   %
   % Call the MEX library.
   %
   try
      [elts] = mice('oscelt_c', state, et, mu );
   catch
      rethrow(lasterror)
   end


