%-Abstract
%
%   CSPICE_SCE2S converts an epoch specified as ephemeris seconds
%   past J2000 (ET) value describing a date to a character string
%   representation of a spacecraft clock value (SCLK).
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
%      sc   the scalar integer NAIF ID of the spacecraft clock whose
%           encoded SCLK value at the epoch 'et' is desired
%
%      et   the double precision scalar or double precision 1xN array of
%           epochs specified as ephemeris seconds past J2000
%
%   the call:
%
%      sclkch = cspice_sce2s( sc, et )
%
%   returns:
%
%      sclkch   the scalar string or NXM character array representation(s)
%               of spacecraft 'sc' clock count  that corresponds to 'et'
%
%               'sclkch' returns with the same vectorization measure (N)
%               as 'et'.
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
%      SC         = -32;
%      SCLK       = '/kernels/voyager2/sclk/vg200004.tsc';
%      event_time = '1979 JUL 05 21:50:21.23379';
%
%      %
%      % Load the SCLK kernel.
%      %
%      cspice_furnsh( SCLK )
%
%      %
%      % Convert the time string to ephemeris time.
%      %
%      et = cspice_str2et( event_time );
%
%      %
%      % Convert the ephemeris time to the corresponding
%      % SCLK string appropriate for this spacecraft
%      %
%      sclkch = cspice_sce2s( SC, et );
%
%      disp( 'Scalar' )
%      txt = sprintf( 'Ephemeris time : %20.8f',  et );
%      disp( txt )
%
%      txt = sprintf( 'SCLK string    : %s', sclkch );
%      disp( txt )
%
%      disp(' ')
%
%      %
%      % Vectorized use, a vector of UTC times.
%      %
%      event_time =  strvcat( '1979 JUL 05 22:50:21.23379', ...
%                             '1979 JUL 05 23:50:21.23379', ...
%                             '1979 JUL 06 00:50:21.23379' );
%
%      %
%      % Convert the time strings to ET.
%      %
%      et = cspice_str2et( event_time );
%
%      %
%      % Convert the ephemeris time to the corresponding
%      % SCLK string appropriate for this spacecraft
%      %
%      sclkch = cspice_sce2s( SC, et );
%
%      disp( 'Vector:' )
%      for i=1:3
%
%         txt = sprintf( 'Ephemeris time : %20.8f',  et(i) );
%         disp( txt )
%
%         txt = sprintf( 'SCLK string    : %s', sclkch(i,:) );
%         disp( txt )
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
%      Scalar
%      Ephemeris time :  -646668528.58223021
%      SCLK string    : 2/20538:39:768
%
%      Vector:
%      Ephemeris time :  -646664928.58223140
%      SCLK string    : 2/20539:54:768
%
%      Ephemeris time :  -646661328.58223248
%      SCLK string    : 2/20541:09:768
%
%      Ephemeris time :  -646657728.58223367
%      SCLK string    : 2/20542:24:768
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine sce2s_c.
%
%   MICE.REQ
%   SCLK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 18-APR-2006, EDW (JPL)
%
%-Index_Entries
%
%   ephemeris time to spacecraft_clock string
%
%-&

function [sclkch] = cspice_sce2s(sc, et)

   switch nargin
      case 2

         sc = zzmice_int(sc);
         et = zzmice_dp(et);

      otherwise
         error ( 'Usage: [_`sclkch`_] = cspice_sce2s(sc, _et_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [sclkch] = mice('sce2s_c',sc, et);
   catch
      rethrow(lasterror)
   end



