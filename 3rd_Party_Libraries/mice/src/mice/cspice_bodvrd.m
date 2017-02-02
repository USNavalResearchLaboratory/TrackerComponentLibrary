%-Abstract
%
%   CSPICE_BODVRD fetches from the kernel pool the double
%   precision values of an item associated with a body.
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
%      bodynm   the scalar string name of the body for which 'item'
%               is requested. 'bodynm' is case-insensitive, and leading
%               and trailing  blanks in 'bodynm' are not significant.
%               Optionally, you may supply the integer ID code for the
%               object as an integer string.  For example both 'MOON'
%               and '301' are legitimate strings that indicate the
%               moon is the body of interest.
%
%      item     the scalar string item name to return. Together, the NAIF
%               ID code of the body and the item name combine to form a
%               kernel variable name, e.g.,
%
%                    'BODY599_RADII'
%                    'BODY401_POLE_RA'
%
%               The values associated with the kernel variable having
%               the name constructed as shown are sought.  Below
%               we'll take the shortcut of calling this kernel variable
%               the "requested kernel variable."
%
%               Note that 'item' *is* case-sensitive.  This attribute
%               is inherited from the case-sensitivity of kernel
%               variable names.
%
%      maxn     the scalar integer defining the maximum number of values
%               the call returns.
%
%   the call:
%
%      values = cspice_bodvrd(body, item, maxn)
%
%   returns:
%
%     values   an array of at most 'maxn' double precision values
%              associated with the requested kernel variable.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      %  Load a set of kernels: an SPK file, a PCK
%      %  file and a leapseconds file. Use a meta
%      %  kernel for convenience.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % When the kernel variable
%      %
%      %    BODY399_RADII
%      %
%      % is present in the kernel pool---normally because a PCK
%      % defining this variable has been loaded (as is the case
%      % here)---the call
%      %
%      values1 = cspice_bodvrd( 'EARTH', 'RADII', 3)
%
%      %
%      % returns the dimension and values associated with the
%      % variable "BODY399_RADII".
%      %
%
%      %
%      % The call lacks case sensitivity in the 'bodynm' variable.
%      %
%     values2 = cspice_bodvrd( 'earth', 'RADII', 3)
%
%      %
%      % The 'item' variable possesses case sensitivity.
%      %
%      try
%
%         %
%         % A call with improper case in 'item' will fail.
%         %
%         values3 = cspice_bodvrd( 'EARTH', 'radii', 3)
%
%      catch
%
%         %
%         % Catch the error, return the error string to the user.
%         %
%         disp( 'Expected error signaled:' )
%         disp( ' ' )
%         disp( lasterr )
%
%      end
%
%      %
%      %  It's always good form to unload kernels after use,
%      %  particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      values1 =
%
%         1.0e+03 *
%
%         6.37814000000000
%         6.37814000000000
%         6.35675000000000
%
%      values2 =
%
%         1.0e+03 *
%
%         6.37814000000000
%         6.37814000000000
%         6.35675000000000
%
%      Expected error signaled:
%
%      SPICE(KERNELVARNOTFOUND): [bodvrd_c->BODVRD] The variable
%      BODY399_radii could not be found in the kernel pool.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine bodvrd_c.
%
%   MICE.REQ
%   KERNEL.REQ
%   NAIF_IDS.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   fetch constants for a body from the kernel pool
%   physical constants for a body
%
%-&

function [values] = cspice_bodvrd(body, item, maxn)

   switch nargin
      case 3

         body = zzmice_str(body);
         item = zzmice_str(item);
         maxn = zzmice_int(maxn);

      otherwise

         error ( 'Usage:  [values()] = cspice_bodvrd( `body`, `item`, maxn)' )

   end

   try
      [values] = mice( 'bodvrd_c', body, item, maxn);
   catch
      rethrow(lasterror)
   end


