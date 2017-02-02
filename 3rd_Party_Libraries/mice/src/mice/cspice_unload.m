%-Abstract
%
%   CSPICE_UNLOAD unloads a SPICE kernel file (of any type)
%   from MATLAB.
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
%      file   the string scalar or NXM character array of SPICE
%             kernel file names, 'file' (or any kernel listed in 'file')
%             should be one loaded through the interface cspice_furnsh
%
%   the call:
%
%      cspice_unload( file )
%
%      removes the file and all associated data from the kernel
%      sub-system. If file is a meta-text kernel, the sub-system
%      unloads all files listed in the kernel.
%
%      Note: a cspice_unload call deletes ALL kernel variables except
%      those loaded into the kernel pool via a cspice_furnsh kernel
%      load  call, i.e. cspice_unload erases kernel variables placed
%      in the pool by the pool functions: cspice_pipool, cspice_pdpool,
%      and cspice_pcpool.
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
%      try
%         values = cspice_bodvrd( 'EARTH', 'RADII', 3);
%         disp('Expected result, found kernel data')
%      catch
%         disp('ERROR: Unexpected result, no kernel data found')
%      end
%
%      %
%      %  Now unload the kernel and try again.
%      %
%      cspice_unload( 'standard.tm' )
%
%      try
%         values = cspice_bodvrd( 'EARTH', 'RADII', 3);
%         disp('ERROR: Unexpected result, found kernel data')
%      catch
%         disp('Expected result, no kernel data found')
%      end
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine unload_c.
%
%   MICE.REQ
%   KERNEL.REQ
%   PCK.REQ
%
%-Version
%
%   -Mice Version 1.1.0, 17-DEC-2008, EDW (JPL)
%
%      Added zzmice_str call on input 'file' to convert string cells to
%      character arrays if 'file' has type string cells. Properly
%      identified 'file' as a vectorizable string/character array.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   Unload a SPICE kernel
%
%-&

function cspice_unload(file)

   switch nargin
      case 1

         file = zzmice_str( file );

      otherwise

         error ( 'Usage: cspice_unload(_`file`_)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('unload_c', file);
   catch
      rethrow(lasterror)
   end

