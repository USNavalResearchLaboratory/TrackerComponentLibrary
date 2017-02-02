%-Abstract
%
%   CSPICE_SCT2E converts encoded spacecraft clock (`ticks')
%   to ephemeris seconds past J2000 (ET).
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
%      sc       the NAIF ID of the spacecraft clock, whose encoded
%               clock value is represented by 'sclkdp'.
%
%               [1,1] = size(sc); int32 = class(sc)
%
%                 or
%
%               [1,1] = size(sc); double = class(sc)
%
%      sclkdp   the encoding of a clock time(s) in units of ticks since the
%               spacecraft clock start time.
%
%               [1,n] = size(sclkdp); double = class(sclkdp)
%
%   the call:
%
%      et = cspice_sct2e( sc, sclkdp )
%
%   returns:
%
%      et    the epoch in ephemeris seconds past J2000, that corresponds 
%            to 'sclkdp'.
%
%            'et' returns with the same vectorization measure (N)
%            as 'sclkdp'.
%
%            [1,n] = size(et); double = class(et)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load the leapseconds kernel for time conversion.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Assign values for the spacecraft ID (Voyager2),
%      % SCLK kernel.
%      %
%      SC     = -32;
%      SCLK   = '/kernels/voyager2/sclk/vg200004.tsc';
%      sclkdp = 985327965.0;
%
%      %
%      % Load the SCLK kernel.
%      %
%      cspice_furnsh( SCLK )
%
%      %
%      % Convert 'sclkdp' for spacecraft 'SC' to ephemeris time.
%      %
%      et = cspice_sct2e( SC, sclkdp );
%
%      %
%      % Convert the ephemeris time to a UTC calendar string.
%      %
%      utc = cspice_et2utc( et, 'C', 3 );
%
%      disp( 'Scalar:' )
%      txt = sprintf( 'SCLKDP: %16.6f', sclkdp );
%      disp( txt )
%
%      txt = sprintf( 'ET    : %16.6f', et );
%      disp( txt )
%
%      txt = sprintf( 'UTC   : %s', utc );
%      disp( txt )
%
%      disp(' ')
%
%      %
%      % Convert a vector of SCLK values.
%      %
%      sclkdp = [ 985327950.0, ...
%                 985553550.0, ...
%                 985901583.0, ...
%                 986447183.0, ...
%                 9136032015.0 ];
%
%      %
%      % Convert the 'sclkdp' vector  for spacecraft 'SC' to
%      % ephemeris time.
%      %
%      et = cspice_sct2e( SC, sclkdp );
%
%      %
%      % Convert the ephemeris time vector to a UTC calendar
%      % strings then output.
%      %
%      utc = cspice_et2utc( et, 'C', 3 );
%
%      disp( 'Vector:' )
%      for i=1:5
%         txt = sprintf( 'SCLKDP: %16.6f', sclkdp(i) );
%         disp( txt )
%
%         txt = sprintf( 'ET    : %16.6f', et(i) );
%         disp( txt )
%
%         txt = sprintf( 'UTC   : %s', utc(i,:) );
%         disp( txt )
%
%         disp(' ')
%
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
%      Scalar:
%      SCLKDP: 985327965.000000
%      ET    : -646668527.682229
%      UTC   : 1979 JUL 05 21:50:22.134
%
%      Vector:
%      SCLKDP: 985327950.000000
%      ET    : -646668528.582228
%      UTC   : 1979 JUL 05 21:50:21.234
%
%      SCLKDP: 985553550.000000
%      ET    : -646654992.592098
%      UTC   : 1979 JUL 06 01:35:57.224
%
%      SCLKDP: 985901583.000000
%      ET    : -646634110.627325
%      UTC   : 1979 JUL 06 07:23:59.189
%
%      SCLKDP: 986447183.000000
%      ET    : -646601374.651195
%      UTC   : 1979 JUL 06 16:29:35.165
%
%      SCLKDP: 9136032015.000000
%      ET    : -157626068.501020
%      UTC   : 1995 JAN 03 02:57:50.315
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine sct2e_c.
%
%   MICE.REQ
%   SCLK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 04-SEP-2012, EDW (JPL)
%
%      Edit to call example in I/O to correct form.
%
%   -Mice Version 1.0.0, 18-APR-2006, EDW (JPL)
%
%-Index_Entries
%
%   spacecraft_clock ticks to ephemeris time
%
%-&

function [et] = cspice_sct2e(sc,sclkdp)

   switch nargin
      case 2

         sc     = zzmice_int(sc);
         sclkdp = zzmice_dp(sclkdp);

      otherwise
         error( 'Usage: [_et_] = cspice_sct2e(sc, _sclkdp_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [et] = mice('sct2e_c',sc, sclkdp);
   catch
      rethrow(lasterror)
   end
