%-Abstract
%
%   MICE_DTPOOL returns descriptive data about a kernel pool variable
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
%     name   the scalar string or NxM character array of names of variables
%            whose values are to be returned
%
%   the call:
%
%      value = mice_dtpool(name)
%
%   returns:
%
%      value   the scalar or 1xN array of structures associating a type
%              and size to a pool variable name.  Each structure
%              contains the fields:
%
%              found   returns as true if the variable 'name' exists
%                      in the pool; false if not
%
%              n      a double precision scalar describing number of values
%                     associated with 'name'. If 'name' does not exist in the
%                     pool, 'n' returns with the value 0
%
%              type   a single character indicating the type of the variable
%                     associated with 'name'
%
%                     'C' if the data is character data
%                     'N' if the data is numeric
%                     'X' if there is no variable name in the pool
%
%              'value' returns with the same vectorization
%               measure (N) as 'name'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a leapsecond kernel.
%      %
%      cspice_furnsh('standard.tm' )
%
%      %
%      % Check for the variables defined in the leapseconds kernel
%      % and names probably (hopefully) not in the kernel pool.
%      %
%      lmpoolNames  = strvcat(              ...
%                    'DELTET/DELTA_T_A',    ...
%                    'DELTET/K',            ...
%                    'DELTET/EB',           ...
%                    'DELTET/M',            ...
%                    'ECHO419',             ...
%                    'DELTET/DELTA_AT',     ...
%                    'EVERLASTING_GOBSTOPPER' );
%
%      value = mice_dtpool( lmpoolNames )
%
%      num   = size(value,2);
%
%      %
%      % Return the indices for which 'found' has value true.
%      %
%      ind = find( [value.found]  );
%
%      for i = 1:size(ind,2)
%
%         fprintf( 'Variable name : %s\n'  , lmpoolNames(ind(i),:)  )
%         fprintf( 'Variable size : %d\n'  , value( ind(i) ).n      )
%         fprintf( 'Variable type : %s\n\n', value( ind(i) ).type   )
%
%      end
%
%      %
%      % Use setdiff to calculate the complement between the 'ind' indices
%      % and the [1:num] set.
%      %
%      ind   = setdiff( [1:num], ind );
%
%      for i = 1:size(ind,2)
%         fprintf('Unable to find variable name : %s\n', lmpoolNames(ind(i),:))
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Variable name : DELTET/DELTA_T_A
%      Variable size : 1
%      Variable type : N
%
%      Variable name : DELTET/K
%      Variable size : 1
%      Variable type : N
%
%      Variable name : DELTET/EB
%      Variable size : 1
%      Variable type : N
%
%      Variable name : DELTET/M
%      Variable size : 2
%      Variable type : N
%
%      Variable name : DELTET/DELTA_AT
%      Variable size : 48
%      Variable type : N
%
%      Unable to find variable name : ECHO419
%      Unable to find variable name : EVERLASTING_GOBSTOPPER
%
%-Particulars
%
%   A sister version of this routine exists named cspice_dtpool that returns
%   the structure field data as separate arguments.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dtpool_c.
%
%-Version
%
%   -Mice Version 1.0.0, 07-MAR-2007, EDW (JPL)
%
%-Index_Entries
%
%   return summary information about a kernel pool variable
%
%-&

function [value] = mice_dtpool(name)

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: [_value_] = mice_dtpool(_`name`_)' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [value] = mice('dtpool_s',name);
   catch
      rethrow(lasterror)
   end





