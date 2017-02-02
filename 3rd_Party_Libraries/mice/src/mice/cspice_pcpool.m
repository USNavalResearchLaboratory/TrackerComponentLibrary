%-Abstract
%
%   CSPICE_PCPOOL provides toolkit programmers a method for
%   programmatically inserting character data into the
%   kernel pool.
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
%              supplied in the array 'cvals'. 'name' is restricted to a length
%              of 32 characters or less.
%
%              [1,m] = size(name); char = class(name)
%
%      cvals   values to load into the kernel pool subsystem with the assigned
%              variable name 'name'.
%
%              [n,m] = size(cvals); char = class(cvals)
%
%   the call:
%
%       cspice_pcpool( name, cvals)
%
%   returns:
%
%      Inserts the variable 'name' into the kernel pool with values as
%      defined in 'cvals'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define the parameters for the string array.
%      %
%      PCPOOL_DIM    = 10;
%      PCPOOL_VAR    = 'pcpool_array';
%      PCPOOL_VAL_TMP= 'pcpool_val';
%      START         = 1;
%
%      %
%      % Populate the 'pcpool_arr' array with values. Initialize
%      % a string array with a string of the correct length.
%      % Note: MATLAB requires the property all strings within
%      % that array have the same length.
%      %
%      pcpool_arr = strvcat( 'n_pcpool_val' );
%
%      %
%      % Fill 'pcpool_arr' with PCPOOL_DIM entries of the
%      % form "n_pcpool_val".
%      %
%      for n=0:PCPOOL_DIM-1
%
%         pcpool_arr(n+1,:) = [ num2str(n) '_' PCPOOL_VAL_TMP ];
%
%      end
%
%      %
%      % Insert the array data into the kernel pool
%      % with variable name 'pcpool_array'.
%      %
%      cspice_pcpool( PCPOOL_VAR, pcpool_arr)
%
%      %
%      % Retrieve the variable's associated values in
%      % array 'cvals'.
%      %
%      cvals = cspice_gcpool( PCPOOL_VAR, START, PCPOOL_DIM );
%
%      %
%      % Check we found the expected variable, and ensure
%      % the expected values.
%      %
%      if ( ~isempty(cvals) )
%
%         txt = sprintf( 'Found array variable %s with entries:', PCPOOL_VAR );
%         disp(txt)
%
%         n_elements = size( cvals, 1);
%
%         for n=1:n_elements
%            txt = sprintf( '   %s', cvals(n,:) );
%            disp(txt)
%         end
%
%      else
%
%         txt = sprintf( 'Failed to find %s in the kernel pool', PCPOOL_VAR );
%         disp(txt)
%
%      end
%
%      %
%      % Clear the kernel pool.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Found array variable pcpool_array with entries:
%         0_pcpool_val
%         1_pcpool_val
%         2_pcpool_val
%         3_pcpool_val
%         4_pcpool_val
%         5_pcpool_val
%         6_pcpool_val
%         7_pcpool_val
%         8_pcpool_val
%         9_pcpool_val
%
%-Particulars
%
%   Kernel pool variable names are restricted to a length of 32
%   characters or less.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pcpool_c.
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
%      Added zzmice_str call on inputs 'name' and 'cvals' to convert
%      string cells to character arrays if 'name' or 'cvals' have
%      type string cells. Added proper markers for usage string
%      variable types.
%
%   -Mice Version 1.0.0, 24-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   Set the value of a character kernel pool variable
%
%-&

function cspice_pcpool( name, cvals )

   switch nargin
      case 2

         name  = zzmice_str(name);
         cvals = zzmice_str(cvals);

      otherwise

         error ( 'Usage: cspice_pcpool( `name`, `cvals`)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('pcpool_c', name, cvals );
   catch
      rethrow(lasterror)
   end



