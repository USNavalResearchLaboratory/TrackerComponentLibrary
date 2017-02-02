%-Abstract
%
%   CSPICE_J2100 returns the value for the Julian Date of
%   2100 JAN 01 12:00:00 (2100 JAN 1.5).
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
%      No input required.
%
%   the call:
%
%      j2100 = cspice_j2100
%
%   returns:
%
%      j2100   the value 2488070.0, the Julian Date corresponding to
%              2100 JAN 01 12:00:00 (2100 JAN 1.5).
%
%              [1,1] = size(j2100); double = class(j2100)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      >> j2100 = cspice_j2100
%
%      j2100 =
%
%           2488070
%
%      >> sprintf( 'J2100 epoch: %10.3f', cspice_j2100 )
%
%      ans =
%
%      J2100 epoch: 2488070.000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine j2100_c.
%
%   MICE.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 11-JUN-2013, EDW (JPL)
%
%       I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   julian date of 2100 jan 1.5
%
%-&

function [j2100] = cspice_j2100

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [j2100] = cspice_j21000' )

   end

   %
   % Call the MEX library.
   %
   try
      [j2100] =  mice('j2100_c');
   catch
      rethrow(lasterror)
   end

