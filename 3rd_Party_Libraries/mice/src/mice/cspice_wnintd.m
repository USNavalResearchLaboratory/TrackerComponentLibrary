%-Abstract
%
%   CSPICE_wnintd returns the window array intersection of two double precision
%   window arrays.
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
%      a   SPICE window containing zero or more intervals.
%
%          [2l,1] = size(a); double = class(a)
%
%      b   SPICE window containing zero or more intervals.
%
%          [2m,1] = size(b); double = class(b)
%
%   the call:
%
%      c = cspice_wnintd( a, b )
%
%   returns:
%
%      c   the double precision window intersection (in the SPICE sense)
%          of 'a' and 'b'
%
%          'c' can overwrite 'a' or 'b'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      SPK1 = 'de405_2000-2050.bsp';
%      SPK2 = 'jup100.bsp';
%
%      %
%      % Retrieve the coverage for body 3 from SPK1
%      %
%      cov1 = cspice_spkcov( SPK1, 3, 10 );
%      fprintf( 'cov1 =\n' )
%      fprintf( '   %16.8f\n', cov1)
%
%   MATLAB outputs:
%
%      cov1 =
%          -43135.81608719
%         1577880064.18391323
%
%      %
%      % Retrieve the coverage for body 3 from SPK2
%      %
%      cov2 = cspice_spkcov( SPK2, 3, 10 );
%      fprintf( 'cov2 =\n' )
%      fprintf( '   %16.8f\n', cov2)
%
%   MATLAB outputs:
%
%      cov2 =
%         -825768000.00000000
%         752241600.00000000
%
%      Perform a windows array intersection on 'cov1' and 'cov2'
%
%      cov3 = cspice_wnintd( cov1, cov2 );
%      fprintf( 'cov3 =\n' )
%      fprintf( '   %16.8f\n', cov3)
%
%   MATLAB outputs:
%
%     cov3 =
%         -43135.81608719
%        752241600.00000000
%
%      The output can overwrite the input.
%
%      cov1 = cspice_wnintd( cov1, cov2 );
%      fprintf( 'cov1 =\n' )
%      fprintf( '   %16.8f\n', cov1)
%
%   MATLAB outputs:
%
%     cov1 =
%         -43135.81608719
%        752241600.00000000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine wnintd_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 09-JUL-2007, EDW (JPL)
%
%-Index_Entries
%
%   intersect two d.p. windows
%
%-&

function [c] = cspice_wnintd( a, b )

   switch nargin

      case 2

         a    = zzmice_win(a);
         b    = zzmice_win(b);

      otherwise

         error ( 'Usage: [c] = cspice_wnintd( a, b )' )

   end

%
% Call the windows routine, add to 'a' and 'b' the space needed for
% the control segments.
%
   try
      [c] = mice('wnintd_c', [zeros(6,1); a], [zeros(6,1); b] );
   catch
      rethrow(lasterror)
   end





