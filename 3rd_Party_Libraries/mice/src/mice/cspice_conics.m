%-Abstract
%
%   CSPICE_CONICS determines the state (position, velocity) of an orbiting
%   body from a set of elliptic, hyperbolic, or parabolic orbital elements.
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
%      elts    a double precision 8-vector or 8xN array containing the conic
%              osculating elements describing the orbit of a body around
%              a primary. The elements are, in order:
%
%                 RP      Perifocal distance.
%                 ECC     Eccentricity.
%                 INC     Inclination.
%                 LNODE   Longitude of the ascending node.
%                 ARGP    Argument of periapse.
%                 M0      Mean anomaly at epoch.
%                 T0      Epoch.
%                 MU      Gravitational parameter.
%
%                 Units are km, rad, rad/sec, km**3/sec**2.
%
%                 The epoch T0 is given in ephemeris seconds past J2000.
%                 T0 is the instant at which the state of the body is
%                 specified by the elements.
%
%      et      the double precision scalar or 1XN-vector of ephemeris
%              time(s) at which to determine the state of the orbiting body
%
%   the call:
%
%      state = cspice_conics(elts, et)
%
%   returns
%
%      state   a double precision Cartesian 6-vector or 6xN array
%              representing the state (position and velocity) of
%              the body at time 'et' in kilometers and kilometers-per-second
%              (the first three components of 'state' represent the x-,
%              y-, and z-components of the body's position; the last three
%              components form the corresponding velocity vector)
%
%              'state' returns with the same vectorization measure (N) as
%              'elts' and 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example (1):
%
%      %
%      % Calculate the perturbation between the
%      % state elements of the Moon at some time as determined
%      % from SPK data and the corresponding state elements
%      % determined from propagation of osculating elements.
%      %
%      % Load the meta kernel listing the needed SPK, PCK, LSK
%      % kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Convert the time of interest, provided as a string, to ephemeris time
%      %
%      et = cspice_str2et( 'Dec 25, 2007' );
%
%      %
%      % Call cspice_spkezr to retrieve the Moon state
%      % w.r.t. the earth in the 'J2000' frame. Use 'NONE' as aberration
%      % correction since we are comparing geometric states.
%      %
%      [state, ltime] = cspice_spkezr( 'Moon', et, 'J2000', 'NONE', 'EARTH' );
%
%      %
%      % cspice_oscelt requires body mass information, so load a
%      % mass PCK kernel that contains gravitation constants.
%      %
%      cspice_furnsh( '/kernels/gen/pck/masses3.tpc' )
%
%      %
%      % Read the gravitational parameter for Earth.
%      %
%      mu = cspice_bodvrd( 'EARTH', 'GM', 1 );
%
%      %
%      % Execute the cspice_oscelt call to convert the state 6-vector
%      % to the osculating elements 8-vector, 'elts', at 'et'. The osculating
%      % elements are relative to the same frame as 'state'.
%      %
%      % The elements describe the nominal orbit the Moon would follow
%      % in the absence of all other bodies in the solar system and
%      % and all non-gravitational forces.
%      %
%      % Note: the cspice_bodvrd call returns data as arrays, so
%      % to access the gravitational parameter (the only value in
%      % the array), we use 'mu(1)'.
%      %
%      elts = cspice_oscelt( state, et, mu(1) );
%
%      %
%      % Now, select a time one week from the initial epoch.
%      %
%      later = et + 7. * cspice_spd;
%
%      %
%      % Use the osculating elements to calculate the state vector
%      % of the Moon at the 'later' epoch.
%      %
%      later_state = cspice_conics( elts, later );
%
%      %
%      % Now retrieve the Moon's state at time 'later' from SPK
%      % data.
%      %
%      [state, ltime] = cspice_spkezr('Moon', later, 'J2000', 'NONE', 'EARTH');
%
%      %
%      % Display the absolute diff between the vector output by
%      % cspice_conics and the state vector returned by cspice_spkezr.
%      %
%      pert = later_state - state;
%
%      txt = sprintf( 'Perturbation in     x: %16.8f', pert(1) );
%      disp( txt )
%
%      txt = sprintf( 'Perturbation in     y: %16.8f', pert(2) );
%      disp( txt )
%
%      txt = sprintf( 'Perturbation in     z: %16.8f', pert(3) );
%      disp( txt )
%
%      txt = sprintf( 'Perturbation in dx/dt: %16.8f', pert(4) );
%      disp( txt )
%
%      txt = sprintf( 'Perturbation in dy/dt: %16.8f', pert(5) );
%      disp( txt )
%
%      txt = sprintf( 'Perturbation in dz/dt: %16.8f', pert(6) );
%      disp( txt )
%
%   MATLAB outputs:
%
%      Perturbation in     x:   -7488.81617036
%      Perturbation in     y:     397.60470311
%      Perturbation in     z:     195.74584983
%      Perturbation in dx/dt:      -0.03615259
%      Perturbation in dy/dt:      -0.00127924
%      Perturbation in dz/dt:      -0.00201456
%
%   Example (2):
%
%      %
%      % Calculate the magnitude of the perturbation between the
%      % position and velocity vectors of the Moon w.r.t. earth as
%      % calculated from cspice_conics and as retrieved from an SPK file.
%      %
%      % Load the meta kernel listing the needed SPK, PCK, LSK
%      % kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Convert the time of interest, provided as a string, to ephemeris time
%      %
%      et1 = cspice_str2et( 'Jan 1, 2007' );
%
%      %
%      % Make the cspice_spkezr call to retrieve the state of the
%      % Moon w.r.t. the earth in J2000. Use 'NONE' as aberration
%      % correction since we are comparing geometric states.
%      %
%      [state1, ltime] = cspice_spkezr( 'Moon', et1, 'J2000', 'NONE', 'EARTH' );
%
%      %
%      % cspice_oscelt requires body mass information, so load a
%      % mass PCK kernel that contains gravitation constants.
%      %
%      cspice_furnsh( '/kernels/gen/pck/masses3.tpc' )
%
%      %
%      % Read the gravitational parameter for Earth.
%      %
%      mu    = cspice_bodvrd( 'EARTH', 'GM', 1 );
%
%      elts1 = cspice_oscelt( state1, et1, mu(1) );
%
%      %
%      % We want to propagate the osculating elements in 'elts1'
%      % by N time steps. Create an array of N copies of 'elts1'
%      %
%      N     = 30;
%      elts  = repmat( elts1, 1, N );
%
%      %
%      % Create an array of N ephemeris times in steps of one day (measured
%      % in seconds) from 'et1'.
%      %
%      et             = [1:N]*cspice_spd + et1;
%
%      twobody        = cspice_conics( elts, et );
%      [state, ltime] = cspice_spkezr( 'Moon', et, 'J2000', 'NONE', 'EARTH' );
%      utc            = cspice_et2utc( et, 'C', 0 );
%
%      for n=1:N
%         txt = sprintf(                                       ...
%                '%s perturbation: ||r|| %10.4f, ||v|| %6.4f', ...
%                 utc(n,:)                                   , ...
%                 norm( state(1:3,n) - twobody(1:3,n) )      , ...
%                 norm( state(4:6,n) - twobody(4:6,n) )            );
%         disp( txt )
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
%                       ... partial output ...
%
%      2007 JAN 12 00:00:00 perturbation: ||r||  5011.4764, ||v|| 0.0282
%      2007 JAN 13 00:00:00 perturbation: ||r||  7828.6919, ||v|| 0.0381
%      2007 JAN 14 00:00:00 perturbation: ||r|| 11573.2356, ||v|| 0.0498
%      2007 JAN 15 00:00:00 perturbation: ||r|| 16336.4334, ||v|| 0.0628
%      2007 JAN 16 00:00:00 perturbation: ||r|| 22123.4631, ||v|| 0.0765
%      2007 JAN 17 00:00:00 perturbation: ||r|| 28830.2006, ||v|| 0.0902
%      2007 JAN 18 00:00:00 perturbation: ||r|| 36232.8928, ||v|| 0.1033
%      2007 JAN 19 00:00:00 perturbation: ||r|| 43994.5246, ||v|| 0.1154
%
%                                ...
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine conics_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   state from conic elements
%
%-&

function [state] = cspice_conics( elts, et )

   switch nargin
      case 2

         elts = zzmice_dp(elts);
         et   = zzmice_dp(et);

      otherwise

         error ( 'Usage: [_state(6)_] = cspice_conics( _elts(8)_, _et_ )' )

   end

   %
   % Call the MEX library.
   %
   try
      [state] = mice('conics_c', elts, et);
   catch
      rethrow(lasterror)
   end



