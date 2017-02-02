%-Abstract
%
%   CSPICE_SCENCD encodes a character representation of spacecraft
%   clock time to the corresponding double precision number.
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
%               whose time is being encoded
%
%      sclkch   the scalar string or NXM character array representation
%               of spacecraft 'sc' clock count
%
%   the call:
%
%      sclkdp = cspice_scencd( sc, sclkch )
%
%   returns:
%
%      sclkdp   the double precision scalar or double precision 1xN array
%               encoding(s) of 'sclkch' for 'sc'
%
%               'sclkdp' returns with the same vectorization measure (N)
%                as 'sclkch'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Assign values for the spacecraft ID (Voyager2),
%      % the SCLK kernel, and a double precision
%      % encodings of SCLK strings
%      %
%      SC     = -32;
%      SCLK   = '/kernels/voyager2/sclk/vg200004.tsc';
%      timein = 985327950.0;
%
%      %
%      % Load the kernel files.
%      %
%      cspice_furnsh( SCLK )
%
%      %
%      % Convert the Voyager encoded SCLK to an
%      % SCLK string.
%      %
%      sclkch = cspice_scdecd( SC, timein );
%
%      %
%      % Convert the SCLK string to double precision form.
%      % The output value should match the original.
%      %
%      sclkdp = cspice_scencd( SC, sclkch );
%
%      disp( 'Scalar:' )
%
%      txt = sprintf( 'Original: %20.8f', timein );
%      disp( txt )
%
%      txt = sprintf( ['SCLKCH  : ' sclkch] );
%      disp( txt )
%
%      txt = sprintf( 'Decoded : %20.8f', sclkdp );
%      disp( txt )
%
%      disp( ' ' )
%
%      %
%      % Convert a vector of SCLK values.
%      %
%      timein = [ 985327950.0, ...
%                 985553550.0, ...
%                 985901583.0, ...
%                 986447183.0, ...
%                 9136032015.0 ];
%
%      %
%      % Convert the SCLK double precision values to the string
%      % representation, then convert to the dp form. As before, the
%      % output value should match the original.
%      %
%      sclkch = cspice_scdecd( SC, timein );
%      sclkdp = cspice_scencd( SC, sclkch );
%
%      disp( 'Vector:' )
%      for i=1:5
%
%         txt = sprintf( 'Original: %20.8f', timein(i) );
%         disp( txt )
%
%         txt = sprintf( ['SCLKCH  : ' sclkch(i,:) ] );
%         disp( txt )
%
%         txt = sprintf( 'Decoded : %20.8f', sclkdp(i) );
%         disp( txt )
%
%         disp( ' ' )
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
%      Original:   985327950.00000000
%      SCLKCH  : 2/20538:39:768
%      Decoded :   985327950.00000000
%
%      Vector:
%      Original:   985327950.00000000
%      SCLKCH  : 2/20538:39:768
%      Decoded :   985327950.00000000
%
%      Original:   985553550.00000000
%      SCLKCH  : 2/20543:21:768
%      Decoded :   985553550.00000000
%
%      Original:   985901583.00000000
%      SCLKCH  : 2/20550:37:001
%      Decoded :   985901583.00000000
%
%      Original:   986447183.00000000
%      SCLKCH  : 2/20561:59:001
%      Decoded :   986447183.00000000
%
%      Original:  9136032015.00000000
%      SCLKCH  : 5/04563:00:001
%      Decoded :  9136032015.00000000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine scencd_c.
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
%   encode spacecraft_clock
%
%-&

function [sclkdp] = cspice_scencd(sc, sclkch)

   switch nargin
      case 2

         sc     = zzmice_int(sc);
         sclkch = zzmice_str(sclkch);

      otherwise
         error ( 'Usage: [_sclkdp_] = cspice_scencd(sc, _`sclkch`_)' )
   end

   %
   % Call the MEX library.
   %
   try
      [sclkdp] = mice('scencd_c',sc, sclkch);
   catch
      rethrow(lasterror)
   end





