%-Abstract
%
%   CSPICE_PDPOOL inserts double precision data into the kernel pool.
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
%      name    name of the kernel pool variable to associate with the values
%              supplied in the array 'dvals'. 'name' is restricted to a length
%              of 32 characters or less.
%
%              [1,m] = size(name); char = class(name)
%
%      dvals   values to load into the kernel pool sub-system with the assigned
%              variable name 'name'.
%
%              [n,1] = size(dvals); double = class(dvals)
%
%   the call:
%
%       cspice_pdpool( name, dvals)
%
%   returns:
%
%      Inserts the variable 'name' into the kernel pool with values as
%      defined in 'dvals'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define the parameters for the double array.
%      %
%      PDPOOL_DIM  = 9;
%      PDPOOL_VAR  = 'pdpool_array';
%      START       = 1;
%
%      %
%      % Populate the 'pdpool_arr' array with PDPOOL_DIM values.
%      %
%      pdpool_arr = [0:PDPOOL_DIM-1]';
%
%      %
%      % Insert the array data into the kernel pool
%      % with variable name 'pipool_array'.
%      %
%      cspice_pdpool( PDPOOL_VAR, pdpool_arr)
%
%      %
%      % Retrieve the variable's associated values in
%      % array 'dvals'.
%      %
%      dvals = cspice_gdpool( PDPOOL_VAR, START, PDPOOL_DIM );
%
%      %
%      % Check we found the expected variable, and ensure
%      % the expected values.
%      %
%      if ( ~isempty(dvals) )
%
%         txt = sprintf( 'Found array variable %s with entries:', PDPOOL_VAR );
%         disp(txt)
%
%         n_elements = size( dvals );
%
%         txt = sprintf( '   %f\n', dvals );
%         disp(txt)
%
%      else
%
%         txt = sprintf( 'Failed to find %s in the kernel pool',  PDPOOL_VAR  );
%         disp(txt)
%
%      end
%
%      %
%      % Clear the kernel pool.
%      %
%      cspice_clpool
%
%   MATLAB outputs:
%
%      Found array variable pdpool_array with entries:
%         0.000000
%         1.000000
%         2.000000
%         3.000000
%         4.000000
%         5.000000
%         6.000000
%         7.000000
%         8.000000
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pdpool_c.
%
%   MICE.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.1.1, 12-MAR-2012, EDW (JPL), SCK (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%      Added mention of the length restriction on the kernel pool variable
%      name 'name'.
%
%   -Mice Version 1.1.0, 23-FEB-2009, EDW (JPL)
%
%      Added zzmice_str call on input 'name' to convert string cells to
%      character arrays if 'name' has type string cells. Added proper
%      markers for usage string variable types.
%
%   -Mice Version 1.0.0, 24-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   Set the value of a d.p. kernel pool variable
%
%-&

function cspice_pdpool( name, dvals )

   switch nargin
      case 2

         name  = zzmice_str(name);
         dvals = zzmice_dp(dvals);

      otherwise

         error ( 'Usage: cspice_pdpool( `name`, dvals)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('pdpool_c', name, dvals );
   catch
      rethrow(lasterror)
   end



