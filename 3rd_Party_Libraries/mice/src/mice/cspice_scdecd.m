%-Abstract
%
%   CSPICE_SCDECD converts a double precision encoding of spacecraft
%   clock time into a string representation.
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
%      sc       the integer scalar NAIF ID of the spacecraft clock
%               whose encoded clock value is represented by
%               'sclkdp'
%
%      sclkdp   the double precision scalar or double precision 1xN array
%               encoding of a clock time(s) in units of ticks since the
%               spacecraft clock start time
%
%   the call:
%
%      sclkch = cspice_scdecd( sc, sclkdp )
%
%   returns:
%
%      sclkch   the scalar string or NXM character array representation(s)
%               of the clock count 'sclkdp' for 'sc'
%
%               'sclkch' returns with the same vectorization measure (N)
%                as 'sclkdp'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Assign values for the spacecraft ID (Voyager2),
%      % SCLK kernels, and an SCLK time string.
%      %
%      SC     = -32;
%      SCLK   = '/kernels/voyager2/sclk/vg200004.tsc';
%      sclkin =  '2/20538:39:768';
%
%      %
%      % Load the SCLK kernel.
%      %
%      cspice_furnsh( SCLK )
%
%      %
%      % Convert the Voyager SCLK strings to the
%      % corresponding double precision value.
%      %
%      timein = cspice_scencd( SC, sclkin );
%
%      %
%      % Convert the double precision value of
%      % the SCLK count back to string. The output
%      % string should nearly match the original
%      % with regards to roundoff and minus any
%      % embedded spaces.
%      %
%      sclkch = cspice_scdecd( SC, timein );
%
%      disp( 'Scalar:' )
%      txt = sprintf( 'Original: %s',  sclkin );
%      disp( txt )
%
%      txt = sprintf( 'Encoded : %20.8f',  timein );
%      disp( txt )
%
%      txt = sprintf( 'Decoded : %s', sclkch );
%      disp( txt )
%
%      %
%      % Convert a vector of SCLK strings. Define a set of strings.
%      %
%      sclkin =  strvcat( '2/20538:39:768' , ...
%                         '2/20543:21:768' , ...
%                         '2/20550:37'     , ...
%                         '2/20561:59'     , ...
%                         '5/04563:00:001'  );
%
%      %
%      % Convert the SCLK strings to the dp representation,
%      % then convert to the string form. As before, the
%      % output value should nearly match the original.
%      %
%      timein = cspice_scencd( SC, sclkin );
%      sclkch = cspice_scdecd( SC, timein );
%
%      disp(' ')
%
%      disp( 'Vector:' )
%      for i=1:5
%
%         txt = sprintf( 'Original: %s',  sclkin(i,:) );
%         disp( txt )
%
%         txt = sprintf( 'Encoded : %20.8f',  timein(i) );
%         disp( txt )
%
%         txt = sprintf( 'Decoded : %s', sclkch(i,:) );
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
%      Scalar:
%      Original: 2/20538:39:768
%      Encoded :   985327950.00000000
%      Decoded : 2/20538:39:768
%
%      Vector:
%      Original: 2/20538:39:768
%      Encoded :   985327950.00000000
%      Decoded : 2/20538:39:768
%
%      Original: 2/20543:21:768
%      Encoded :   985553550.00000000
%      Decoded : 2/20543:21:768
%
%      Original: 2/20550:37
%      Encoded :   985901583.00000000
%      Decoded : 2/20550:37:001
%
%      Original: 2/20561:59
%      Encoded :   986447183.00000000
%      Decoded : 2/20561:59:001
%
%      Original: 5/04563:00:001
%      Encoded :  9136032015.00000000
%      Decoded : 5/04563:00:001
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine scdecd_c.
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
%   decode spacecraft_clock
%
%-&

function [sclkch] = cspice_scdecd(sc, sclkdp)

   switch nargin
      case 2

        sc     = zzmice_int(sc);
        sclkdp = zzmice_dp(sclkdp);

      otherwise

         error ( 'Usage: [_`sclkch`_] = cspice_scdecd(sc, _sclkdp_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [sclkch] = mice('scdecd_c',sc,sclkdp);
   catch
      rethrow(lasterror)
   end



