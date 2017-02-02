%-Abstract
%
%   CSPICE_KCLEAR clears the KEEPER system: unload all kernels, clears
%   the kernel pool, and re-initialize the system.
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
%   The call:
%
%      cspice_kclear
%
%      Re-initialize the KEEPER system.
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
%   This routine allows you re-initialize the KEEPER system with
%   a single call.  The KEEPER system is the kernel management system
%   underlying the set of Mice APIs
%
%      cspice_furnsh
%      cspice_ktotal
%      cspice_kdata
%      cspice_kinfo
%      cspice_kclear
%      cspice_unload
%
%   This routine unloads all kernels from their kernel-type-specific
%   kernel management subsystems (SPKBSR, CKBSR, etc.), clears the
%   kernel pool, clears KEEPER's internal file database, and re-sets
%   the watch status for the kernel variables used to load kernels
%   via meta-kernels.
%
%   This capability, though implemented in Fortran, is particularly
%   relevant to SPICE implementations such as Mice, for which the
%   state of the KEEPER system persists after any Mice-based MATLAB
%   script is run. Successive runs of Mice-based scripts may perform
%   in unexpected ways when scripts access data loaded during runs of
%   previous scripts.
%
%   Cleaning up after such programs using explicit unload_c commands is
%   tedious and error-prone.  One call to this routine sets the
%   KEEPER system to its initial state, preventing unintentional
%   interaction between scripts via KEEPER's state.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine kclear_c.
%
%   MICE.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 01-DEC-2006, EDW (JPL)
%
%-Index_Entries
%
%   Re-initialize the keeper system
%   Clear the keeper system
%   Unload all kernels
%
%-&

function cspice_kclear

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: cspice_kclear' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('kclear_c');
   catch
      rethrow(lasterror)
   end



