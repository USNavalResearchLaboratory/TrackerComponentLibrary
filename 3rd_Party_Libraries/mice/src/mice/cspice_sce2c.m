%-Abstract
%
%   CSPICE_SCE2C converts ephemeris seconds past J2000 (ET) to
%   continuous encoded spacecraft clock ("ticks").  Non-integral
%   tick values may be returned.
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
%      et   the scalar or N-vector of double precision epochs,
%           specified as ephemeris seconds past J2000
%
%   the call:
%
%      sclkdp = cspice_sce2c( sc, et )
%
%   returns:
%
%      sclkdp   the double precision scalar or double precision 1xN array
%               of encoded SCLK value(s) corresponding to 'et' for 'sc'
%
%               'sclkdp' returns with the same vectorization
%                measure (N) as 'et'.
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
%      SC         = -32;
%      SCLK       = '/kernels/voyager2/sclk/vg200004.tsc';
%      event_time = '1979 JUL 05 21:50:21.23379';
%
%      %
%      % Load the SCLK file.
%      %
%      cspice_furnsh( SCLK )
%
%      %
%      % Convert the time string to ephemeris time.
%      %
%      et = cspice_str2et( event_time );
%
%      %
%      % Convert the ephemeris time to the encoded SCLK
%      % format.
%      %
%      sclkdp = cspice_sce2c( SC, et );
%      txt    = sprintf( ' %16.6f', sclkdp );
%      disp( txt )
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
%      % Convert the 'et' array to the encoded
%      % spacecraft clock.
%      %
%      sclkdp = cspice_sce2c( SC, et );
%
%      for i=1:3
%         txt = sprintf( ' %16.6f', sclkdp(i) );
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
%      Scalar:
%       985327949.999971
%
%      Vector:
%       985387950.043701
%       985447950.087433
%       985507950.131163
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine sce2c_c.
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
%   ephemeris time to continuous spacecraft_clock ticks
%
%-&

function [sclkdp] = cspice_sce2c(sc, et)

   switch nargin
      case 2

         sc = zzmice_int(sc);
         et = zzmice_dp(et);

      otherwise
         error ( 'Usage: [_sclkdp_] = cspice_sce2c(sc, _et_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [sclkdp] = mice('sce2c_c',sc, et);
   catch
      rethrow(lasterror)
   end



