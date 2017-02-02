%-Abstract
%
%   CSPICE_KTOTAL returns the current number of kernels loaded
%   via the KEEPER interface that are of a specified type.
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
%     kind   a scalar string list of types of kernels to count when
%            checking the loaded kernels. 'kind' should consist of a list of
%            words of kernels to examine.  Recognized types are
%
%                 SPK  --- All SPK files are counted in the total.
%                 CK   --- All CK files are counted in the total.
%                 PCK  --- All binary PCK files are counted in the
%                          total.
%                 EK   --- All EK files are counted in the total.
%                 TEXT --- All text kernels that are not meta-text.
%                          kernels are included in the total.
%                 META --- All meta-text kernels are counted in the
%                          total.
%                 ALL  --- Every type of kernel is counted in the
%                          total.
%
%            'kind' lacks case sensitivity. The cspice_ktotal algorithm ignores
%            words in 'kind' if not one of those listed above.
%
%            See the Examples section for illustrations of the
%            use of kind.
%
%   the call:
%
%      count = cspice_ktotal( kind )
%
%   returns:
%
%      count   a double precision scalar describing the number of kernels
%              loaded through cspice_furnsh belonging to the list
%              specified by 'kind'
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%     %
%     % Load the standard meta kernel, retrieve the number of
%     % loaded kernels.
%     %
%     cspice_furnsh( 'standard.tm' )
%
%     n   = cspice_ktotal( 'ALL' );
%     txt = sprintf('Count of loaded kernels before cspice_kclear call: %d', n);
%     disp( txt )
%
%   MATLAB outputs:
%
%     Count of loaded kernels before cspice_kclear call: 4
%
%   The expected result counting standard.tm and the three kernels
%   named in the meta kernel.
%
%     %
%     % Clear the KEEPER system, retrieve the number of loaded
%     % after the clear.
%     %
%     cspice_kclear
%
%     n   = cspice_ktotal( 'ALL' );
%     txt = sprintf('Count of loaded kernels after cspice_kclear call: %d', n);
%     disp( txt )
%
%   MATLAB outputs:
%
%     Count of loaded kernels after cspice_kclear call: 0
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ktotal_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 06-MAY-2009, EDW (JPL)
%
%      Added MICE.REQ reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 01-DEC-2006, EDW (JPL)
%
%-Index_Entries
%
%   Number of loaded kernels of a given type
%
%-&

function [count] = cspice_ktotal( kind )

   switch nargin
      case 1

         kind = zzmice_str(kind);

      otherwise

         error ( 'Usage: [count] = cspice_ktotal(`kind`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [count] = mice( 'ktotal_c', kind );

      %
      % Convert the integers returned from the interface to double precision
      % in case a user includes the return arguments in a calculation
      % with other doubles.
      %
      count = zzmice_dp(count);

   catch
      rethrow(lasterror)
   end



