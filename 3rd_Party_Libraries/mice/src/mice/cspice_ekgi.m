%-Abstract
%
%   CSPICE_EKGI returns an element of integer data from a
%   specified row in a specified column of the set of rows matching
%   the previous cspice_ekfind SELECT query.
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
%      selidx   the scalar integer index for a column of interest
%               satisfying the SELECT clause, the column indices
%               range from 1 to number of columns in the SELECT clause.
%
%      row      the scalar integer index for a row in the column
%               identified by 'selidx', the column indices
%               range from 1 to 'nmrows' where 'nmrows' equals the total
%               number of rows satisfying the SELECT clause.
%
%      elment   the scalar integer index for an element of
%               the data at the 'selidx','row' position; a scalar
%               value at 'selidx', 'row' has 'elment' value one.
%
%   the call:
%
%      [ idata, null, found] = cspice_ekgi( selidx, row, elment )
%
%   returns:
%
%      idata    the integer value of the requested element at
%               data location 'selidx', 'row', 'elment'.
%
%      null     a scalar boolean indicating if 'idata' has a null value.
%
%      found    a scalar boolean indicating whether the specified
%               value at 'selidx', 'row', 'elment' was found.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Assign an EK file to load..
%      %
%      EK = 'test_file.ek';
%
%      %
%      % Load the EK.
%      %
%      cspice_furnsh( EK )
%
%      %
%      % Assume the file test_file.ek contains the table 'scalar_2',
%      % and that "scalar_2' has the column named 'i_col_1' of integer
%      % values.
%      %
%      % Define a set of constraints to perform a query on all
%      % loaded EK files (the SELECT clause). In this case select
%      % the column "i_col_1" from table "scalar_2."
%      %
%      query = 'Select i_col_1 from scalar_2 order by row_no';
%
%      %
%      % Query the EK system for data rows matching the
%      % SELECT constraints.
%      %
%      [ nmrows, ok, errmsg ] = cspice_ekfind( query );
%
%      %
%      % Check whether an error occurred while processing the
%      % SELECT clause. If so, output the error message.
%      %
%      if ( ok )
%         printf( 'SELECT clause error: %s\n', errmsg );
%      end
%
%      %
%      % Loop over each row found matching the query.
%      %
%      for rowno = 1:nmrows
%
%         %
%         % Fetch the integer data. We know the query returned
%         % one column and the column contains only scalar data,
%         % so the index of all elements is 1.
%         %
%         selidx = 1;
%         eltidx = 1;
%
%         %
%         % Use cspice_ekgi to retrieve the value from
%         % row/column position.
%         %
%         [ idata, isnull, found ] = cspice_ekgi( selidx, ...
%                                                 rowno,  ...
%                                                 eltidx );
%
%         %
%         % Output the value, if non-null data exist at the
%         % requested position.
%         %
%         if  ~isnull
%            fprintf( 'Integer data: %d\n', idata );
%         end
%
%      end
%
%      %
%      % Clear the kernel pool and database. Note, you don't normally
%      % unload an EK after a query, rather at the end of a program.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Integer data: 2000001
%      Integer data: 2000002
%      Integer data: 2000003
%      Integer data: 2000004
%      Integer data: 2000005
%      Integer data: 2000006
%      Integer data: 2000007
%      Integer data: 2000008
%      Integer data: 2000009
%      Integer data: 2000010
%      Integer data: 2000011
%      Integer data: 2000012
%      Integer data: 2000013
%      Integer data: 2000014
%      Integer data: 2000015
%      Integer data: 2000016
%      Integer data: 2000017
%      Integer data: 2000018
%      Integer data: 2000019
%      Integer data: 2000020
%
%-Particulars
%
%   Suppose a SELECT clause return data consisting of three columns (N=3)
%   and four rows (M=4):
%
%              col 1    col 2    col 3
%
%      row 1   val_11   val_12   val_13
%      row 2   val_21   val_22   val_23
%      row 3   val_31   val_32   val_33
%      row 4   val_41   val_42   val_43
%
%   with "col 2" and "col 3" containing scalar integer data and "val_42"
%   containing a vector of K integers.
%
%   Retrieving the data elements depends on the values for the index set
%   "selidx," "row," and "elment."
%
%   Use the set
%
%      'selidx' = 2, 'row' = 3, 'elment' = 1
%
%   to fetch scalar "val_32."
%
%   Use the set
%
%      'selidx' = 3, 'row' = 4, 'elment' = 1
%
%   to fetch scalar "val_43."
%
%   Use the set
%
%      'selidx' = 2, 'row' = 4, 'elment' = K
%
%   to fetch the final element of vector "val_42"
%
%   'elment' is allowed to exceed the number of elements in the column
%   entry; if it does, 'found' returns as false.  This allows the caller
%   to read data from the column entry in a loop without checking the
%   number of available elements first.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ekgi_c.
%
%   MICE.REQ
%   EK.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 10-MAY-2011, EDW (JPL)
%
%      "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.0, 10-APR-2010, EDW (JPL)
%
%-Index_Entries
%
%   fetch element from integer column entry
%
%-&

function [ idata, null, found] = cspice_ekgi( selidx, row, elment )

   switch nargin
      case 3

         selidx = zzmice_int(selidx);
         row    = zzmice_int(row);
         elment = zzmice_int(elment);

      otherwise

         error ( [ 'Usage: [ idata, null, found] = ' ...
                   'cspice_ekgi( selidx, row, elment )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [ idata, null, found] = mice('ekgi_c', selidx, row, elment );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      null  = zzmice_logical(null);
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end


