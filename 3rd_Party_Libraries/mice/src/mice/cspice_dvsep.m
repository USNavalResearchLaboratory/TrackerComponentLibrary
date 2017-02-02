%-Abstract
%
%   CSPICE_DVSEP calculates the time derivative of the separation angle
%   between states.
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
%      s1   defining a SPICE state(s);
%
%           [6,n] = size(s1); double = class(s1)
%
%              s1 = (r1, dr1 ).
%                         --
%                         dt
%
%      s2   defining a second SPICE state(s);
%
%           [6,n] = size(s2); double = class(s2)
%
%              s2 = (r2, dr2 ).
%                        --
%                        dt
%
%      An implicit assumption exists that 's1' and 's2' are specified
%      in the same reference frame. If this is not the case, the numerical
%      result has no meaning.
%
%   the call:
%
%      dvsep = cspice_dvsep( s1, s2)
%
%   returns:
%
%      dvsep   time derivative(s) of the angular separation between 's1' and
%              's2'.
%
%              [1,n] = size(dvsep); double = class(dvsep)
%
%              'dvsep' returns with the same measure of vectorization (N)
%              as 's1' and 's2'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load SPK, PCK, and LSK kernels, use a meta kernel for convenience.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % An arbitrary time.
%      %
%      BEGSTR = 'JAN 1 2009';
%      et     = cspice_str2et( BEGSTR );
%
%      %
%      % Calculate the state vectors sun to Moon, sun to earth at ET.
%      %
%      [statee, ltime] = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN' );
%      [statem, ltime] = cspice_spkezr('MOON',  et, 'J2000', 'NONE', 'SUN' );
%
%      %
%      % Calculate the time derivative of the angular separation of
%      % the earth and Moon as seen from the sun at ET.
%      %
%      dsept = cspice_dvsep( statee, statem );
%      fprintf( 'Time derivative of angular separation, rads/sec: %.10e\n', ...
%                                                                   dsept )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Time derivative of angular separation, rads/sec: 3.8121193603e-09
%
%-Particulars
%
%   In this discussion, the notation
%
%      < V1, V2 >
%
%   indicates the dot product of vectors V1 and V2. The notation
%
%      V1 x V2
%
%   indicates the cross product of vectors V1 and V2.
%
%   To start out, note that we need consider only unit vectors,
%   since the angular separation of any two non-zero vectors
%   equals the angular separation of the corresponding unit vectors.
%   Call these vectors U1 and U2; let their velocities be V1 and V2.
%
%   For unit vectors having angular separation
%
%      THETA
%
%   the identity
%
%      || U1 x U1 || = ||U1|| * ||U2|| * sin(THETA)                (1)
%
%   reduces to
%
%      || U1 x U2 || = sin(THETA)                                  (2)
%
%   and the identity
%
%      | < U1, U2 > | = || U1 || * || U2 || * cos(THETA)           (3)
%
%   reduces to
%
%      | < U1, U2 > | = cos(THETA)                                 (4)
%
%   Since THETA is an angular separation, THETA is in the range
%
%      0 : Pi
%
%   Then letting s be +1 if cos(THETA) > 0 and -1 if cos(THETA) < 0,
%     we have for any value of THETA other than 0 or Pi
%
%
%                                2          1/2
%      cos(THETA) = s * ( 1 - sin (THETA)  )                       (5)
%
%   or
%
%                                2          1/2
%      < U1, U2 > = s * ( 1 - sin (THETA)  )                       (6)
%
%
%   At this point, for any value of THETA other than 0 or Pi,
%   we can differentiate both sides with respect to time (T)
%   to obtain
%
%                                                    2        -1/2
%      < U1, V2 > + < V1, U2 > =    s * (1/2)(1 - sin (THETA))
%
%                                 * (-2) sin(THETA)*cos(THETA)
%
%                                 * d(THETA)/dT                   (7a)
%
%
%   Using equation (5), and noting that s = 1/s, we can cancel
%   the cosine terms on the right hand side
%
%                                                    -1
%      < U1, V2 > + < V1, U2 > =    (1/2)(cos(THETA))
%
%                                 * (-2) sin(THETA)*cos(THETA)
%
%                                 * d(THETA)/dT                   (7b)
%
%   With (7b) reducing to
%
%      < U1, V2 > + < V1, U2 > = - sin(THETA) * d(THETA)/dT        (8)
%
%   Using equation (2) and switching sides, we obtain
%
%      || U1 x U2 || * d(THETA)/dT  =  - < U1, V2 > - < V1, U2 >   (9)
%
%   or, provided U1 and U2 are linearly independent,
%
%      d(THETA)/dT = ( - < U1, V2 > - < V1, U2 > ) / ||U1 x U2||  (10)
%
%   Note for times when U1 and U2 have angular separation 0 or Pi
%   radians, the derivative of angular separation with respect to
%   time doesn't exist. (Consider the graph of angular separation
%   with respect to time; typically the graph is roughly v-shaped at
%   the singular points.)
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dvsep_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   time derivative of angular separation
%
%-&

function [dvsep] = cspice_dvsep(s1, s2)

   switch nargin
      case 2

         s1 = zzmice_dp(s1);
         s2 = zzmice_dp(s2);

      otherwise

         error ( 'Usage: [_dvsep_] = cspice_dvsep(_s1(6)_, _s2(6)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [dvsep] = mice('dvsep_c', s1, s2);
   catch
      rethrow(lasterror)
   end



