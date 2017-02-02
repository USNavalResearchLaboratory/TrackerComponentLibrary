%-Abstract
%
%   CSPICE_ET2UTC converts an input time from ephemeris seconds
%   past J2000 to Calendar, Day-of-Year, or Julian Date format, UTC.
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
%      et       the double precision scalar or 1xN array of ephemeris
%               time expressed as ephemeris seconds past J2000
%
%      format   the scalar string format flag describing the output time
%               string, it may be any of the following:
%
%                  'C'      Calendar format, UTC
%
%                  'D'      Day-of-Year format, UTC
%
%                  'J'      Julian Date format, UTC
%
%                  'ISOC'   ISO Calendar format, UTC
%
%                  'ISOD'   ISO Day-of-Year format, UTC
%
%      prec     the scalar integer number of decimal places of precision to
%               which fractional seconds (for Calendar and Day-of-Year
%               formats) or days (for Julian Date format) are to be
%               computed
%
%   the call:
%
%      utcstr = cspice_et2utc( et, format, prec )
%
%   returns:
%
%      utcstr    the scalar string or NXM character array of output time
%                strings equivalent to the input epoch 'et', in the specified
%                'format'
%
%                'utcstr' returns with the same vectorization measure (N)
%                as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define an arbitrary ephemeris time.
%      %
%      et     = -527644192.5403653;
%      format = 'J';
%      prec   = 6;
%      SIZE   = 5;
%
%      %
%      % Load a leapseconds kernel.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Convert the ephemeris time to Julian Date
%      % 'format'. Define precision to 6 decimal
%      % places.
%      %
%      utcstr = cspice_et2utc( et, format, prec );
%      disp( 'Scalar:' )
%
%      txt = sprintf( 'ET              : %12.4f', et );
%      disp(txt)
%
%      txt = sprintf( 'Converted output: %s', utcstr );
%      disp( txt )
%
%      %
%      % Create an array of ephemeris times beginning
%      % at -527644192.5403653 with graduations of 10000.0
%      % ephemeris seconds.
%      %
%      et     = [0:(SIZE-1)]*10000. -527644192.5403653;
%      format = 'C';
%
%      %
%      % Convert the array of ephemeris times 'et' to an
%      % array of UTC strings, 'utcstr', in calendar
%      % 'format'.
%      %
%      utcstr= cspice_et2utc( et, format, prec );
%
%      disp( ' ' )
%      disp( 'Vector:' )
%
%      for n=1:SIZE
%
%         txt = sprintf( 'ET              : %12.4f', et(n) );
%         disp( txt )
%
%         txt = sprintf( 'Converted output: %s', utcstr(n,:) );
%         disp( txt )
%
%         disp(' ' )
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
%      ET              : -527644192.5404
%      Converted output: JD 2445438.006415
%
%      Vector:
%      ET              : -527644192.5404
%      Converted output: 1983 APR 13 12:09:14.274000
%
%      ET              : -527634192.5404
%      Converted output: 1983 APR 13 14:55:54.274001
%
%      ET              : -527624192.5404
%      Converted output: 1983 APR 13 17:42:34.274001
%
%      ET              : -527614192.5404
%      Converted output: 1983 APR 13 20:29:14.274002
%
%      ET              : -527604192.5404
%      Converted output: 1983 APR 13 23:15:54.274002
%
%-Particulars
%
%   Use of this routine requires a loaded leapseconds kernel.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine et2utc_c.
%
%   MICE.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   ephemeris time to utc
%
%-&

function [utcstr] = cspice_et2utc(et, format, prec )

   switch nargin
      case 3

         et     = zzmice_dp(et);
         format = zzmice_str(format);
         prec   = zzmice_int(prec);

      otherwise

         error ( 'Usage: [_`utcstr`_] = cspice_et2utc(_et_, `format`, prec)' )

   end

   %
   % Call the MEX library.
   %
   try
      [utcstr] = mice('et2utc_c',et,format,prec);
   catch
      rethrow(lasterror)
   end


