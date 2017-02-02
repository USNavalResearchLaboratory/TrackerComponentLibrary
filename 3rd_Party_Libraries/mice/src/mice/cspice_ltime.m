%-Abstract
%
%   CSPICE_LTIME computes the transmit (or receive) time of
%   a signal at a specified target, given the receive (or transmit)
%   time at a specified observer. The elapsed time between transmit
%   and receive is also returned.
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
%      etobs   a double precision scalar or 1xN vector defining the
%              epoch in ephemeris seconds (TDB) of a signal at some
%              observer
%
%      obs     the scalar integer NAIF ID code of the observer
%
%      dir     a character string pictograph defining the
%              direction the signal travels, to target from
%              observer "->", or from the target to the
%              observer "<-"
%
%      targ    the scalar integer NAIF ID code of the target
%
%   the call:
%
%      [ettarg, elapsd] = cspice_ltime(etobs, obs, dir, targ)
%
%   returns:
%
%      ettarg   the double precision scalar or 1XN vector defining the
%               epoch at which the electromagnetic signal is "at" the
%               target body, expressed in ephemeris seconds (TDB)
%
%                  Note 'ettarg' is computed using only Newtonian
%                  assumptions about the propagation of light.
%
%      elapsd   the double precision scalar or 1XN vector defining the
%               measure of ephemeris seconds (TDB) between transmission
%               and receipt of the signal
%
%                  elapsd = abs( etobs - ettarg )
%
%               'ettarg' and 'elapsd' return with the same
%               vectorization measure (N) as 'etobs'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      %  Load an SPK, PCK, and leapseconds kernel
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Suppose a signal originates from Earth towards the
%      % the Jupiter system barycenter. Define the NAIF IDs
%      % for the observer, Earth (399), the target, Jupiter
%      % barycenter (5), and time of interest.
%      %
%      OBS      = 399;
%      TARGET   = 5;
%      TIME_STR = 'July 4, 2004';
%
%      %
%      %  Convert the transmission time to ET.
%      %
%      et = cspice_str2et( TIME_STR);
%
%      %
%      %  Determine the arrival time and the time for propagation.
%      %
%      [arrive, ltime] = cspice_ltime( et, OBS, '->', TARGET);
%
%      %
%      %  Convert the arrival time (ET) to UTC.
%      %
%      arrive_utc = cspice_et2utc( arrive, 'C', 3 );
%
%      %
%      %  Output the results.
%      %
%      txt = sprintf( 'Transmission at (UTC)       : %s', TIME_STR );
%      disp(txt)
%
%      txt = sprintf( 'The signal arrived at (UTC) : %s', arrive_utc );
%      disp(txt)
%
%      txt = sprintf( 'Time for propagation (secs) : %16.4f', ltime );
%      disp(txt)
%      disp( ' ' )
%
%      %
%      % Now assume the signal originated at Jupiter barycenter,
%      % received by Earth at TIME_STR. Determine the transmission
%      % time and the time for propagation.
%      %
%      [receive, ltime] = cspice_ltime( et, OBS, '<-', TARGET);
%
%      %
%      % Convert the reception time (ET) to UTC.
%      %
%      receive_utc = cspice_et2utc( receive, 'C', 3 );
%
%      %
%      %  Output the results.
%      %
%      txt = sprintf( 'Reception at (UTC)          : %s', TIME_STR );
%      disp(txt)
%
%      txt = sprintf( 'The signal sent at (UTC)    : %s', receive_utc );
%      disp(txt)
%
%      txt = sprintf( 'Time for propagation (secs) : %16.4f', ltime );
%      disp(txt)
%
%   MATLAB outputs:
%
%      Transmission at (UTC)       : July 4, 2004
%      The signal arrived at (UTC) : 2004 JUL 04 00:48:38.717
%      Time for propagation (secs) :        2918.7170
%
%      Reception at (UTC)          : July 4, 2004
%      The signal sent at (UTC)    : 2004 JUL 03 23:11:21.248
%      Time for propagation (secs) :        2918.7524
%
%-Particulars
%
%     None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ltime_c.
%
%   MICE.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   Compute uplink and downlink light time
%
%-&

function [ettarg, elapsd] = cspice_ltime(etobs, obs, dir, targ)

   switch nargin
      case 4

         etobs = zzmice_dp(etobs);
         obs   = zzmice_int(obs);
         targ  = zzmice_int(targ);
         dir   = zzmice_str(dir);

      otherwise

         error ( ['Usage: [_ettarg_, _elapsd_] = ' ...
                  'cspice_ltime( _etobs_, obs, `dir`, targ)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [ettarg, elapsd] = mice('ltime_c',etobs, obs, dir, targ);
   catch
      rethrow(lasterror)
   end

