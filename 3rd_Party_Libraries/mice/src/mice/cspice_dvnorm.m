
%-Abstract
%
%   CSPICE_DVNORM returns the derivative of the vector norm of a 3-vector.
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
%      state   6-vector(s), the second three components of the vector(s)
%              being the derivatives of the first three with respect to
%              some scalar.
%
%                 state =  ( x, dx )
%                               --
%                               ds
%
%              A common form for 'state' would contain position and
%              velocity.
%
%              [6,n] = size(state); double = class(state)
%
%   the call:
%
%      dvnorm = cspice_dvnorm(state)
%
%   returns:
%
%      dvnorm   the value(s) of d||x|| corresponding to 'state'.
%                               ------
%                               ds
%
%                                    1/2         2    2    2  1/2
%               Where ||x|| = < x, x >    =  ( x1 + x2 + x3 )
%
%
%                         v = ( dx1, dx2, dx3 )
%                               ---  ---  ---
%                               ds   ds   ds
%
%                     d||x||   < x, v >
%                    ------ =   ------    =  < xhat, v >
%                      ds             1/2
%                              < x, x >
%
%             'dvnorm' returns with the same measure of vectorization (N)
%             as 'state'.
%
%              [1,n] = size(dvnorm); double = class(dvnorm)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Create several 6-vectors (6x1 arrays) with the structure
%      %
%      %   y = |  x  |
%      %       |     |
%      %       |  dx |
%      %       |  -- |
%      %       |  ds |
%      %
%      % where 'x' is a 3-vector (3x1 array).
%      %
%
%      %
%      % Create 'y' with 'x' of varying magnitudes. Use 'x'
%      % and '-x' to define the derivative as parallel and
%      % anti-parallel.
%      %
%      mag = [ -4, 4, 12 ];
%
%      x   = [ 1, sqrt(2), sqrt(3 ) ]';
%
%      y   = [ [x * 10^mag(1);  x], ...
%              [x * 10^mag(2); -x], ...
%              [  zeros(3,1);  x * 10^mag(3) ] ];
%
%      %
%      % Calculate the derivative of the vector norms with respect
%      % to 's'.
%      %
%      dvnorm = cspice_dvnorm( y );
%
%      fprintf( 'Parallel x, dx/ds         : %f\n', dvnorm(1) )
%      fprintf( 'Anti-parallel x, dx/ds    : %f\n', dvnorm(2) )
%      fprintf( 'Zero vector x, large dx/ds: %f\n', dvnorm(3) )
%
%   Matlab outputs:
%
%      Parallel x, dx/ds         : 2.449490
%      Anti-parallel x, dx/ds    : -2.449490
%      Zero vector x, large dx/ds: 0.000000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dvnorm_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 09-NOV-2012 (EDW)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 10-MAY-2010 (EDW)
%
%-Index_Entries
%
%   derivative of 3-vector norm
%
%-&

function [dvnorm] = cspice_dvnorm(state)

   switch nargin
      case 1

         state = zzmice_dp(state);

      otherwise

         error ( 'Usage: [_dvnorm_] = cspice_dvnorm(_state(6)_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [dvnorm] = mice('dvnorm_c', state);
   catch
      rethrow(lasterror)
   end
