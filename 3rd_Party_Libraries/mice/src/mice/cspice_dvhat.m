%-Abstract
%
%   CSPICE_DVHAT calculates the unit vector corresponding to a state or states
%   and the derivative of the unit vector.
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
%      s1   a double precision 6x1-array or 6xN array defining a
%           state or states;
%
%              s1 = (r1, dr1 ).
%                         --
%                         dt
%
%   the call:
%
%      dvhat = cspice_dvhat(s1)
%
%   returns:
%
%      dvhat   a double precision 6x1 array or 6xN array containing the unit
%              vector(s) pointing in the direction of the position component(s)
%              of 's1' and the derivative of the unit vector with respect
%              to time;
%
%              dvhat = [u, du ] where u =   r1
%                          --             -----
%                          dt             ||r1||
%
%              'dvhat' returns with the same measure of vectorization (N)
%              as 's1'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Suppose that 'state' gives the apparent state of a body with
%   respect to an observer.  This routine can be used to compute the
%   instantaneous angular rate of the object across the sky as seen
%   from the observers vantage.
%
%      %
%      % Load SPK, PCK, and LSK kernels, use a meta kernel for convenience.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define an arbitrary epoch, convert the epoch to ephemeris time.
%      %
%      EPOCH = 'Jan 1 2009';
%      et    = cspice_str2et( EPOCH );
%
%      %
%      % Calculate the state of the moon with respect to the earth-moon
%      % barycenter in J2000, corrected for light time and stellar aberration
%      % at 'et'.
%      %
%      target   = 'MOON';
%      frame    = 'J2000';
%      abcorr   = 'LT+S';
%      observer = 'EARTH BARYCENTER';
%
%      [ state, ltime ] = cspice_spkezr( target, et, frame, abcorr, observer );
%
%      %
%      % Calculate the unit vector of 'state' and the derivative of the
%      % unit vector.
%      %
%      ustate = cspice_dvhat( state )
%
%      %
%      % Calculate the instantaneous angular velocity from the magnitude of the
%      % derivative of the unit vector.
%      %
%      %   v = r x omega
%      %
%      %   ||omega|| = ||v||  for  r . v = 0
%      %               -----
%      %               ||r||
%      %
%      %   ||omega|| = ||v||  for  ||r|| = 1
%      %
%      omega = cspice_vnorm( ustate(4:6) );
%
%      fprintf( 'Instantaneous angular velocity %2.10e rad/sec.\n', omega )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Instantaneous angular velocity 2.4810665797e-06  rad/sec.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dvhat_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 04-MAY-2010, EDW (JPL)
%
%-Index_Entries
%
%   state of a unit vector parallel to a state vector
%
%-&

function [dvhat] = cspice_dvhat(s1)

   switch nargin
      case 1

         s1 = zzmice_dp(s1);

      otherwise

         error ( 'Usage: [_dvhat(6)_] = cspice_dvhat(_s1(6)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [dvhat] = mice('dvhat_c',s1);
   catch
      rethrow(lasterror)
   end
