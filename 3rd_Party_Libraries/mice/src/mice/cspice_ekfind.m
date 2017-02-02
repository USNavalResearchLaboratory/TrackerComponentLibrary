%-Abstract
%
%   CSPICE_EKFIND finds E-kernel data that satisfy a set of constraints.
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
%      query   a string scalar specifying the data to locate from data
%              available in all loaded EK files. The general form of a
%              query general form:
%
%                 SELECT   <column list>
%                 FROM     <table list>
%                 [WHERE    <constraint list>]
%                 [ORDER BY <ORDER BY column list>]
%
%              (WHERE and ORDER BY are optional parameters)
%
%   the call:
%
%      [ nmrows, ok, errmsg] = cspice_ekfind( query )
%
%   returns:
%
%      nmrows   a scalar integer containing the number of rows matching
%               the query.
%
%      ok       a scalar boolean indicating whether the query parsed
%               correctly (TRUE) or not (FALSE).
%
%      errmsg   a string scalar containing a description of the parse
%               error should one occur, otherwise the string returns
%               as blank.
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
%      % Assume the file test_file.ek contains the table 'scalar_2',
%      % and that "scalar_2' contains columns named:
%      %
%      %   c_col_1, d_col_1, i_col_1, t_col_1
%      %
%      % Define a set of constraints to perform a query on all
%      % loaded EK files (the SELECT clause).
%      %
%      query = [ 'Select c_col_1, d_col_1, i_col_1, t_col_1 from ' ...
%                'scalar_2 order by row_no' ];
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
%      % If no error occurred, 'nmrows' contains the number of rows matching
%      % the constraints specified in the query string.
%      %
%      fprintf( 'Number of matching row: %d\n', nmrows )
%
%      %
%      % Clear the kernel pool and database. Note, you don't normally
%      % unload an EK after a query, rather at the end of a program.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Number of matching row: 20
%
%   Load at least one EK kernel prior to calling cspice_ekfind, otherwise
%   an error signals.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ekfind_c.
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
%   find EK data
%   issue EK query
%
%-&

function [ nmrows, ok, errmsg] = cspice_ekfind(query)

   switch nargin
      case 1

         sample = zzmice_str(query);

      otherwise

         error ( [ 'Usage: [ nmrows, ok, `errmsg`] = ' ...
                   'cspice_ekfind( `query` )' ])

   end

   %
   % Call the MEX library.
   %
   try
      [ nmrows, ok, errmsg] = mice('ekfind_c', query );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      ok = zzmice_logical(ok);
   catch
      rethrow(lasterror)
   end


