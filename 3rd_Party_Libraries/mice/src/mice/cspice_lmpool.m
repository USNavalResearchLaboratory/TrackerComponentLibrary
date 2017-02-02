%-Abstract
%
%   CSPICE_LMPOOL loads the variables contained in a text buffer
%   into the kernel pool.
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
%      cvals   a scalar string (1xM array), NxM array of strings, or cell
%              of N strings defining SPICE kernel variable assignments
%              that could serve as a SPICE text kernel.
%
%   the call:
%
%       cspice_lmpool( cvals)
%
%   inserts the variable assignments defined by 'cvals' into the
%   kernel pool subsystem. Once inserted, the user can access the
%   variables using the cspice_gcpool, cspice_gipool, or cspice_gdpool
%   calls.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      lmpoolNames  = {                   ...
%                    'DELTET/DELTA_T_A',  ...
%                    'DELTET/K',          ...
%                    'DELTET/EB',         ...
%                    'DELTET/M',          ...
%                    'DELTET/DELTA_AT'    ...
%                     };
%
%      %
%      % Create a kernel in a text buffer.
%      %
%      textbuf = {                                             ...
%             'DELTET/DELTA_T_A = 32.184',                     ...
%             'DELTET/K         = 1.657D-3',                   ...
%             'DELTET/EB        = 1.671D-2',                   ...
%             'DELTET/M         = ( 6.239996 1.99096871D-7 )', ...
%             'DELTET/DELTA_AT  = ( 10, @1972-JAN-1',          ...
%             '                     11, @1972-JUL-1',          ...
%             '                     12, @1973-JAN-1',          ...
%             '                     13, @1974-JAN-1',          ...
%             '                     14, @1975-JAN-1',          ...
%             '                     15, @1976-JAN-1',          ...
%             '                     16, @1977-JAN-1',          ...
%             '                     17, @1978-JAN-1',          ...
%             '                     18, @1979-JAN-1',          ...
%             '                     19, @1980-JAN-1',          ...
%             '                     20, @1981-JUL-1',          ...
%             '                     21, @1982-JUL-1',          ...
%             '                     22, @1983-JUL-1',          ...
%             '                     23, @1985-JUL-1',          ...
%             '                     24, @1988-JAN-1',          ...
%             '                     25, @1990-JAN-1',          ...
%             '                     26, @1991-JAN-1',          ...
%             '                     27, @1992-JUL-1',          ...
%             '                     28, @1993-JUL-1',          ...
%             '                     29, @1994-JUL-1',          ...
%             '                     30, @1996-JAN-1',          ...
%             '                     31, @1997-JUL-1',          ...
%             '                     32, @1999-JAN-1 )'         ...
%                };
%
%      %
%      % Load the kernel data into the kernel pool.
%      %
%      cspice_lmpool( textbuf )
%
%      %
%      % Ensure the loaded data exists in the kernel pool.
%      % Query the pool for each expected name, size of the
%      % variable with that name, and the type of data
%      % for that name.
%      %
%
%     [found, n, type] = cspice_dtpool( lmpoolNames );
%
%      for i = 1:numel(lmpoolNames)
%
%         if ( found(i) )
%
%            fprintf( ['Found %s, with %i values assigned' ...
%                      ' of data type %s.\n\n'],   ...
%                      char(lmpoolNames(i)), n(i), type(i) )
%
%         end
%
%      end
%
%      %
%      %  It's always good form to unload kernels after use,
%      %  particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Found DELTET/DELTA_T_A, with 1 values assigned of data type N.
%
%      Found DELTET/K, with 1 values assigned of data type N.
%
%      Found DELTET/EB, with 1 values assigned of data type N.
%
%      Found DELTET/M, with 2 values assigned of data type N.
%
%      Found DELTET/DELTA_AT, with 46 values assigned of data type N.
%
%-Particulars
%
%   This routine allows you to store a text kernel in an internal
%   array of your program and load this array into the kernel pool
%   without first storing its contents as a text kernel.
%
%   Kernel pool variable names are restricted to a length of 32
%   characters or less.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine lmpool_c.
%
%   MICE.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 10-FEB-2010, EDW (JPL)
%
%      Added mention of the length restriction on kernel pool variable
%      names.
%
%   -Mice Version 1.0.0, 23-FEB-2009, EDW (JPL)
%
%-Index_Entries
%
%   Load the kernel pool from an internal text buffer
%
%-&

function cspice_lmpool( cvals )

   switch nargin
      case 1

         cvals = zzmice_str( cvals);

      otherwise

         error ( 'Usage: cspice_lmpool( _`cvals`_ )' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('lmpool_c', cvals );
   catch
      rethrow(lasterror)
   end



