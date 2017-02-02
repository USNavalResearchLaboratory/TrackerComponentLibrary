%-Abstract
%
%   CSPICE_EKNELT returns the number of elements in a specified column entry
%   in the current row.
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
%      selidx   the scalar integer index of the column in the
%               SELECT from which to retrieve data. The range of
%               'selidx' is 1:nsel inclusive, where 'nsel' is the
%               number of items in the SELECT clause of the current
%               query.
%
%      row      the scalar integer index of the row containing the element.
%               This number refers to a member of the set of rows
%               matching a query. 'row' must be in the range
%
%                  1:nmrows
%
%               where 'nmrows' is the matching row count returned
%               by cspice_ekfind.
%
%   the call:
%
%       nelt = cspice_eknelt( selidx, row )
%
%   returns:
%
%      nelt    the scalar integer number of elements in the column entry
%              belonging to the specified column in the specified row.
%
%      Null entries in variable-size columns are considered to have size 1.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Assign an EK file to load.
%      %
%      EK = 'test_file.ek';
%
%      %
%      % Load the EK.
%      %
%      cspice_furnsh( EK )
%
%      %
%      % The file "test_file.ek" contains the table 'vector_1', and
%      % 'vector_1' has the column named 'd_col_1', a vector of double
%      % precision values.
%      %
%
%      %
%      % Define a set of constraints to perform a query on all
%      % loaded EK files (the SELECT clause). In this case select
%      % the column "d_col_1" from table "vector_1."
%      %
%      query = 'Select d_col_1 from vector_1 order by row_no';
%
%      %
%      % Query the EK system for data rows matching the
%      % SELECT restraints.
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
%         % Fetch the double precision data. We know the query returned
%         % one column, determine the number of values in the row of
%         % interest.
%         %
%         selidx = 1;
%         nelt   = cspice_eknelt( selidx, rowno);
%
%         %
%         % Use cspice_ekgd to retrieve the value from
%         % row/column position.
%         %
%
%         for eltidx = 1:nelt
%
%            [ ddata, isnull, found ] = cspice_ekgd( selidx, ...
%                                                    rowno,  ...
%                                                    eltidx );
%
%            %
%            % Output the value, if non-null data exist at the
%            % requested position.
%            %
%            if  ~isnull
%               fprintf( 'Double precision data (%d,%d,%d): %f\n', ...
%                        selidx, rowno, eltidx, ddata );
%            end
%
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
%   Matlab outputs:
%
%      Double precision data (1,1,1): 5000101.000000
%      Double precision data (1,1,2): 5000102.000000
%      Double precision data (1,1,3): 5000103.000000
%      Double precision data (1,1,4): 5000104.000000
%      Double precision data (1,2,1): 5000201.000000
%      Double precision data (1,2,2): 5000202.000000
%      Double precision data (1,2,3): 5000203.000000
%      Double precision data (1,2,4): 5000204.000000
%      Double precision data (1,3,1): 5000301.000000
%      Double precision data (1,3,2): 5000302.000000
%      Double precision data (1,3,3): 5000303.000000
%      Double precision data (1,3,4): 5000304.000000
%
%-Particulars
%
%   This routine is meant to be used in conjunction with the EK fetch
%   entry points cspice_ekgc, cspice_ekgd, and cspice_ekgi.  This routine
%   allows the caller of those routines to determine appropriate
%   loop bounds to use to fetch each column entry in the current row.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine eknelt_c.
%
%   MICE.REQ
%   EK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 10-APR-2010, EDW (JPL)
%
%-Index_Entries
%
%   return the number of elements in a column entry
%
%-&

function [nelt] = cspice_eknelt( selidx, row )

   switch nargin

      case 2

         selidx = zzmice_int(selidx);
         row    = zzmice_int(row);

      otherwise

         error ( 'Usage: [nelt] = cspice_eknelt( selidx, row )' )

   end

   try
      [nelt] = mice('eknelt_c', selidx, row );
   catch
      rethrow(lasterror)
   end
