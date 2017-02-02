%-Abstract
%
%   CSPICE_SCS2E converts a spacecraft clock string to ephemeris
%   seconds past J2000 (ET).
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
%      sc       the scalar integer NAIF ID of the spacecraft clock
%               whose clock value is represented by 'sclkch'
%
%      sclkch   the scalar string or NXM character array representation
%               of spacecraft 'sc' clock count ('sclkch' is an absolute
%               spacecraft clock time, so the string should include
%               partition information)
%
%   the call:
%
%      et = cspice_scs2e( sc, sclkch )
%
%   returns:
%
%      et    the double precision scalar or double precision 1xN array
%            in ephemeris seconds past J2000, that corresponds to 'sclkch'
%
%            'et' returns with the same vectorization measure
%             (N) as 'sclkch'.
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
%      % and SCLK kernel.
%      %
%      SC     = -32;
%      SCLK   = '/kernels/voyager2/sclk/vg200004.tsc';
%      sclkch = '2/20538:39:768';
%
%      %
%      % Load the SCLK kernel.
%      %
%      cspice_furnsh( SCLK )
%
%      %
%      % Convert 'sclkch' for spacecraft 'SC' to ephemeris time.
%      %
%      et = cspice_scs2e( SC, sclkch );
%
%      %
%      % Convert the ephemeris time to a UTC calendar string.
%      %
%      utc = cspice_et2utc( et, 'C', 3 );
%
%      disp( 'Scalar:' )
%      txt = sprintf( 'Original:  %s', sclkch );
%      disp( txt )
%
%      txt = sprintf( 'ET      : %20.8f',  et );
%      disp( txt )
%
%      txt = sprintf( 'UTC     : %s', utc );
%      disp( txt )
%
%      disp (' ')
%
%      %
%      % Convert a vector of SCLK strings to ET and
%      % UTC.
%      %
%      sclkch =  strvcat( '2/20538:39:768' , ...
%                         '2/20543:21:768' , ...
%                         '2/20550:37'     , ...
%                         '2/20561:59'     , ...
%                         '5/04563:00:001'  );
%
%      et  = cspice_scs2e( SC, sclkch );
%      utc = cspice_et2utc( et, 'C', 3 );
%
%      disp( 'Vector:' )
%      for i=1:5
%
%         txt = sprintf( 'Original:  %s', sclkch(i,:) );
%         disp( txt )
%
%         txt = sprintf( 'ET      : %20.8f',  et(i) );
%         disp( txt )
%
%         txt = sprintf( 'UTC     : %s', utc(i,:) );
%         disp( txt )
%
%         disp (' ')
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
%      Original:  2/20538:39:768
%      ET      :  -646668528.58222842
%      UTC     : 1979 JUL 05 21:50:21.234
%
%      Vector:
%      Original:  2/20538:39:768
%      ET      :  -646668528.58222842
%      UTC     : 1979 JUL 05 21:50:21.234
%
%      Original:  2/20543:21:768
%      ET      :  -646654992.59209847
%      UTC     : 1979 JUL 06 01:35:57.224
%
%      Original:  2/20550:37
%      ET      :  -646634110.62732494
%      UTC     : 1979 JUL 06 07:23:59.189
%
%      Original:  2/20561:59
%      ET      :  -646601374.65119493
%      UTC     : 1979 JUL 06 16:29:35.165
%
%      Original:  5/04563:00:001
%      ET      :  -157626068.50102001
%      UTC     : 1995 JAN 03 02:57:50.315
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine scs2e_c.
%
%   MICE.REQ
%   SCLK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 18-APR-2006, EDW (JPL)
%
%-Index_Entries
%
%   spacecraft_clock string to ephemeris time
%
%-&

function [et] = cspice_scs2e(sc, sclkch)

   switch nargin
      case 2

         sc     = zzmice_int(sc);
         sclkch = zzmice_str(sclkch);

      otherwise
         error ( 'Usage: [_et_] = cspice_scs2e(sc, _`sclkch`_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [et] = mice('scs2e_c',sc, sclkch);
   catch
      rethrow(lasterror)
   end





