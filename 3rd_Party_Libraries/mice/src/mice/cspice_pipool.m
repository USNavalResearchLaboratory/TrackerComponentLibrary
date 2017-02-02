%-Abstract
%
%   CSPICE_PIPOOL provides toolkit programmers a method for
%   programmatically inserting integer data into the kernel pool.
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
%              supplied in the array 'ivals'. 'name' is restricted to a length
%              of 32 characters or less.
%
%              [1,m] = size(name); char = class(name)
%
%      ivals   values to load into the kernel pool subsystem with the assigned
%              variable name 'name'.
%
%              [n,1] = size(ivals); int32 = class(ivals)
%
%   the call:
%
%       cspice_pipool( name, ivals)
%
%   returns:
%
%      Inserts the variable 'name' into the kernel pool with values as
%      defined in 'ivals'.
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
%      PIPOOL_DIM  = 9;
%      PIPOOL_VAR  = 'pipool_array';
%      START       = 1;
%
%      %
%      % Populate the 'pipool_arr' array with PIPOOL_DIM values.
%      %
%      pipool_arr = int32([0:PIPOOL_DIM-1]');
%
%      %
%      % Insert the array data into the kernel pool
%      % with variable name 'pipool_array'.
%      %
%      cspice_pipool( PIPOOL_VAR, pipool_arr)
%
%      %
%      % Retrieve the variable's associated values in
%      % array 'ivals'.
%      %
%      ivals = cspice_gipool( PIPOOL_VAR, START, PIPOOL_DIM );
%
%      %
%      % Check we found the expected variable, and ensure
%      % the expected values.
%      %
%      if ( ~isempty(ivals) )
%
%         txt = sprintf( 'Found array variable %s with entries:', PIPOOL_VAR );
%         disp(txt)
%
%         n_elements = size( ivals );
%
%         txt = sprintf( '   %d\n', ivals );
%         disp(txt)
%
%      else
%
%         txt = sprintf( 'Failed to find %s in the kernel pool',  PIPOOL_VAR  );
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
%      Found array variable pipool_array with entries:
%         0
%         1
%         2
%         3
%         4
%         5
%         6
%         7
%         8
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pipool_c.
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
%   Set the value of a numeric kernel pool variable
%
%-&

function cspice_pipool( name, ivals )

   switch nargin
      case 2

         name  = zzmice_str(name);
         ivals = zzmice_int(ivals);

      otherwise

         error ( 'Usage: cspice_pipool( `name`, ivals[])' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('pipool_c', name, ivals );
   catch
      rethrow(lasterror)
   end



