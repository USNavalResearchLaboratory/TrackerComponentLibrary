%-Abstract
%
%   CSPICE_CLPOOL clears the kernel pool.
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
%      cspice_clpool
%
%   deletes all variable assignments loaded into the kernel
%   pool.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Place a value into the kernel pool. Recall
%      % the routines for direct insertion
%      % of pool assignments have arrays for input,
%      % but in MATLAB a scalar is a 1x1 array.
%      %
%      cspice_pdpool( 'TEST_VAR', -666. )
%
%      %
%      % Check for the variable assignment to TEST_VAR.
%      % cspice_gdpool returns an empty array if the variable
%      % does not exist in the kernel pool.
%      %
%      dvals = cspice_gdpool( 'TEST_VAR', 0, 1 );
%
%      if ( ~isempty(dvals) )
%         disp( sprintf( 'TEST_VAR value: %f', dvals ) )
%      end
%
%      %
%      % Now clear the kernel pool.
%      %
%      cspice_clpool
%
%      %
%      % Again, check for the TEST_VAR assignment.
%      %
%      dvals = cspice_gdpool( 'TEST_VAR', 0, 1 );
%
%      if ( isempty(dvals)  )
%         disp( 'TEST_VAR not in kernel pool' )
%      end
%
%   MATLAB outputs, after the first cspice_gdpool call:
%
%      TEST_VAR value: -666.000000
%
%   Demonstrating the existence of the assignment on the
%   kernel pool.
%
%   MATLAB outputs, after the second cspice_gdpool call:
%
%      TEST_VAR not in kernel pool
%
%   The variable assignment no longer exists in the kernel pool.
%
%-Particulars
%
%   Note, cspice_clpool deletes ALL pool assignments, including those
%   from cspice_boddef and the cspice_pipool, cspice_pdpool, cspice_pcpool
%   set. Use cspice_unload to remove the assignments loaded from a
%   particular kernel.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine clpool_c.
%
%   MICE.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   CLEAR the pool of kernel variables
%
%-&

function cspice_clpool

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: cspice_clpool' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('clpool_c');
   catch
      rethrow(lasterror)
   end

