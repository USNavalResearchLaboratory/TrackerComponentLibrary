%-Abstract
%
%   CSPICE_CKCOV returns the coverage windows for a specified object
%   in a specified CK file.
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
%      ck        the string name, or cell of strings, of SPICE CK files.
%
%                [1,c] = size(ck), char = class(ck)
%
%                  or
%
%                [1,m] = size(ck), cell = class(ck)
%
%      idcode    the integer CK ID code of an object, normally a
%                spacecraft structure or instrument, for which pointing data
%                are expected to exist in the specified CK file
%
%                [1,1] = size(idcode), int32 = class(idcode)
%
%      needav    a boolean scalar indicating whether to consider segments
%                having angular velocity when determining coverage. When
%                'needav' is true, segments without angular velocity don't
%                contribute to the coverage window; when 'needav' is false,
%                all segments for 'idcode' may contribute to the coverage
%                window.
%
%                [1,1] = size(needav), logical = class(needav)
%
%      level     the string defining the level (granularity) at which
%                the coverage is examined.  Allowed values and corresponding
%                meanings are:
%
%                     'SEGMENT'    The output coverage window contains
%                                  intervals defined by the start and
%                                  stop times of segments for the object
%                                  designated by 'idcode'.
%
%                     'INTERVAL'   The output coverage window contains
%                                  interpolation intervals of segments
%                                  for the object designated by
%                                  'idcode'.  For type 1 segments, which
%                                  don't have interpolation intervals,
%                                  each epoch associated with a pointing
%                                  instance is treated as a singleton
%                                  interval; these intervals are added
%                                  to the coverage window.
%
%                                  All interpolation intervals are
%                                  considered to lie within the segment
%                                  bounds for the purpose of this
%                                  summary:  if an interpolation
%                                  interval extends beyond the segment
%                                  coverage interval, only its
%                                  intersection with the segment
%                                  coverage interval is considered to
%                                  contribute to the total coverage.
%
%                [1,c] = size(level), char = class(level)
%
%      tol       a double precision scalar defining the tolerance value
%                expressed in ticks of the spacecraft clock associated with
%                'idcode'. Before each interval is inserted into the coverage
%                window, the interval is intersected with the segment coverage
%                interval, then and if the intersection is non-empty, it is
%                expanded by 'tol': the left endpoint of the intersection
%                interval is reduced by and the right endpoint is increased by
%                'tol'. Adjusted interval 'tol' endpoints, when expressed as
%                encoded SCLK, never are less than zero ticks. Any intervals
%                that overlap as a result of the expansion are merged.
%
%                The coverage window returned when tol > 0 indicates the
%                coverage provided by the file to the CK readers cspice_ckgpav
%                and cspice_ckgp when that value of 'tol' is passed to them as
%                an input.
%
%                [1,1] = size(tol), double = class(tol)
%
%      timsys    the string scalar indicating the time system used in the
%                output coverage window. 'timsys' may have the
%                values:
%
%                      'SCLK'    Elements of 'cov' are expressed in
%                                encoded SCLK ("ticks"), where the
%                                clock is associated with the object
%                                designated by 'idcode'.
%
%                      'TDB'     Elements of 'cov' are expressed as
%                                seconds past J2000 TDB.
%
%                [1,c] = size(timsys), char = class(timsys)
%
%      room      an integer scalar defining the number of intervals for use
%                as a workspace by the routine. This value should equal at least
%                the number of intervals corresponding to 'idcode' in 'ck'.
%
%                [1,1] = size(room), int32 = class(room)
%
%      cover_i   an optional input describing a either an empty window or a
%                window array created from a previous cspice_ckcov call.
%                Inclusion of this window argument results in an output
%                window consisting of a union of the data retrieved from the
%                'ck' kernels and the data in 'cover_i'.
%
%                [2m,1] = size(cover_i), double = class(cover_i)
%
%                   or
%
%                [0,0] = size(cover_i), double = class(cover_i)
%
%   the call:
%
%      cov = cspice_ckcov(ck, idcode, needav, level, tol, timsys, room, cover_i)
%
%         or
%
%      cov = cspice_ckcov(ck, idcode, needav, level, tol, timsys, room)
%
%   returns:
%
%      cov   a double precision 2Nx1 array containing the coverage for 'idcode'
%            available from 'ck'.  When the coverage level is 'INTERVAL', this
%            is the set of time intervals for which data for 'idcode' are
%            present in the file 'ck', The array 'cov' contains the pairs of
%            endpoints of these intervals.
%
%            Each window defined as a pair of endpoints such that:
%
%               window 1 = cover(1:2)
%               window 2 = cover(3:4)
%               window 3 = cover(5:6)
%                           ...
%               window N = cover(2N-1,2N)
%
%            When the coverage 'level' is 'SEGMENT', 'cov' is computed in a
%            manner similar to that described above, but the coverage intervals
%            used in the computation are those of segments rather than
%            interpolation intervals within segments.
%
%            When 'tol' is > 0, the intervals comprising the coverage window for
%            'idcode' are expanded by 'tol' and any intervals overlapping as a
%            result are merged. The resulting window is returned in 'cov'.
%            The expanded window in no case extends beyond the segment bounds
%            in either direction by more than 'tol'.
%
%            The interval endpoints contained in 'cov' are encoded spacecraft
%            clock times if 'timsys' is 'SCLK'; otherwise the times are
%            converted from encoded spacecraft clock to seconds past J2000 TDB.
%
%            'cov' returns an empty set if 'ck' lacks coverage for 'idcode'. If
%            'cover_i' exists in the argument list, 'cov' returns as a union of
%            the coverage data found in 'ck' and the data in 'cover_i'. 'cov'
%            can overwrite 'cover_i'.
%
%            [2p,1] = size(cov), double = class(cov)
%
%               or
%
%            [0,1] = size(cov), double = class(cov)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Use a simple function to display the CK IDs found in a CK, or set of
%   CKs, and the time coverage of the data corresponding to those IDs.
%   This example calls both cspice_ckobj and cspice_ckcov. In practice,
%   algorithms using cspice_ckobj will also use cspice_ckcov and
%   vice-versa.
%
%   function ckcov_t( CK, SCLK, LEVEL )
%
%         MAXIV  = 100000;
%         WINSIZ = 2 * MAXIV;
%         MAXOBJ = 1000;
%         LSK    = 'naif0010.tls';
%
%         %
%         % Load a leapseconds kernel and the SCLK corresponding to the
%         % input CK.
%         %
%         % Note, neither cspice_ckcov or cspice_ckobj require these
%         % kernels to function. We need these data for output time
%         % conversion.
%         %
%         cspice_furnsh( LSK )
%         cspice_furnsh( SCLK)
%
%         %
%         % Find the set of objects in the CK file.
%         %
%         ids = cspice_ckobj( CK, MAXOBJ );
%
%         %
%         % We want to display the coverage for each object. Loop over
%         % the contents of the ID code set, find the coverage for
%         % each item in the set, and display the coverage.
%         %
%         for i=1:numel(ids)
%
%            %
%            % Extract the coverage data for object 'ids(i)'.
%            %
%            cover    = cspice_ckcov(CK, ids(i), 0, LEVEL, 0.0, 'TDB', WINSIZ);
%            [row,col]= size(cover);
%
%            %
%            % Display a simple banner.
%            %
%            fprintf( '========================================\n')
%            fprintf( 'Coverage for object %d\n', ids(i) )
%
%            %
%            %  'cover' has dimension 2Nx1, where 'row' has the value 2N with
%            %  each window defined as a pair of endpoints such that:
%            %
%            %  window 1 = cover(1:2)
%            %  window 2 = cover(3:4)
%            %  window 3 = cover(5:6)
%            %        ...
%            %  window N = cover(2N-1,2N)
%            %
%            % Loop from 1 to 'row' with stepsize 2.
%            %
%            for j=1:2:row
%
%               %
%               % Convert the endpoints to TDB calendar format time strings
%               % and display them. Pass the endpoints in an array,
%               % so cspice_timout returns an array of time strings.
%               %
%               % Recall a vectorized input has dimension 1xM so transpose
%               % the 'cover' slice.
%               %
%               timstr = cspice_timout( cover(j:j+1)', ...
%                                   'YYYY MON DD HR:MN:SC.### (TDB) ::TDB' );
%               fprintf('Interval: %d\n'  , (j+1)/2 )
%               fprintf('   Start: %s\n'  , timstr(1,:) )
%               fprintf('    Stop: %s\n\n', timstr(2,:) )
%
%            end
%
%         end
%
%         %
%         % Empty the kernel pool.
%         %
%         cspice_kclear
%
%   Example (1):
%
%      >> SCLK = '/kernels/cassini/sclk/cas00101.tsc';
%      >> CK   = { '/kernels/cassini/ck/05357_05362ra.bc', ...
%                  '/kernels/cassini/ck/05362_06002ra.bc'   };
%
%      Output data using the 'INTERVAL' level.
%
%      >> ckcov_t( CK, SCLK, 'INTERVAL' )
%
%   MATLAB outputs:
%
%      ========================================
%      Coverage for object -82000
%      Interval: 1
%         Start: 2005 DEC 23 00:01:07.900 (TDB)
%          Stop: 2005 DEC 23 15:36:55.540 (TDB)
%
%      Interval: 2
%         Start: 2005 DEC 23 15:37:39.539 (TDB)
%          Stop: 2005 DEC 23 16:59:35.508 (TDB)
%
%      Interval: 3
%         Start: 2005 DEC 23 17:00:43.507 (TDB)
%          Stop: 2005 DEC 24 13:55:59.025 (TDB)
%
%      Interval: 4
%         Start: 2005 DEC 24 13:56:19.024 (TDB)
%          Stop: 2005 DEC 24 17:25:42.944 (TDB)
%
%                ... continued ...
%
%      Interval: 24
%         Start: 2005 DEC 31 15:49:11.103 (TDB)
%          Stop: 2006 JAN 01 15:18:34.561 (TDB)
%
%      Interval: 25
%         Start: 2006 JAN 01 15:20:30.560 (TDB)
%          Stop: 2006 JAN 01 16:43:38.528 (TDB)
%
%      Interval: 26
%         Start: 2006 JAN 01 16:45:02.528 (TDB)
%          Stop: 2006 JAN 01 22:52:10.386 (TDB)
%
%      Interval: 27
%         Start: 2006 JAN 01 22:52:38.386 (TDB)
%          Stop: 2006 JAN 02 00:01:02.360 (TDB)
%
%   Example (2):
%
%      Output data using the 'SEGMENT' level.
%
%      >> ckcov_t( CK, SCLK, 'SEGMENT' )
%
%   MATLAB outputs:
%
%      ========================================
%      Coverage for object -82000
%      Interval: 1
%         Start: 2005 DEC 23 00:01:07.900 (TDB)
%          Stop: 2005 DEC 28 00:01:01.130 (TDB)
%
%      Interval: 2
%         Start: 2005 DEC 28 00:01:05.130 (TDB)
%          Stop: 2006 JAN 02 00:01:02.360 (TDB)
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ckcov_c.
%
%   MICE.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 03-APR-2012, EDW (JPL)
%
%      Edits to Example code and comments. No change to Example code
%      functionality.
%
%      Renamed the argument 'size' to 'room'. "size" is a Matlab function
%      name and it's seriously dumb to use a function name word as an argument
%      name.
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%      Explicitly described ID variables as "CK IDs."
%
%   -Mice Version 1.1.0, 29-DEC-2008, EDW (JPL)
%
%      Edited description of 'size'; 'size' now defines the maximum
%      number of intervals for the internal workspace window.
%
%      The 'cover_i' argument may now have the empty array value, [],
%      on input.
%
%      Added range restriction on size.
%
%      Corrected misspellings.
%
%   -Mice Version 1.0.0, 18-JUN-2007, EDW (JPL)
%
%-Index_Entries
%
%   get coverage window for ck object
%
%-&

function [cov] = cspice_ckcov( ck, idcode, needav, level, tol, ...
                                           timsys, room, cover_i )

   switch nargin
      case 7

         ck     = zzmice_str(ck);
         idcode = zzmice_int(idcode);
         needav = zzmice_int(needav);
         level  = zzmice_str(level);
         tol    = zzmice_dp(tol);
         timsys = zzmice_str(timsys);
         room   = zzmice_int(room, [1, int32(inf)/2] );

      case 8

         ck     = zzmice_str(ck);
         idcode = zzmice_int(idcode);
         needav = zzmice_int(needav);
         level  = zzmice_str(level);
         tol    = zzmice_dp(tol);
         timsys = zzmice_str(timsys);
         room   = zzmice_int(room, [1, int32(inf)/2] );
         cover_i= zzmice_win(cover_i);

      otherwise

         error ( ['Usage: [cov] = cspice_ckcov( _`ck`_, idcode, ' ...
                          'needav, level, tol, timsys, room, [cover_i])'] )

   end

%
% The call passed either seven or eight arguments. Branch accordingly.
%
if ( nargin == 7 )

   %
   % Call the MEX library.
   %
   try
      cov = mice('ckcov_c', ck, idcode, needav, level, tol, timsys, room );
   catch
      rethrow(lasterror)
   end


else

   %
   % Call the MEX library.
   %
   try
      cov = mice('ckcov_c', ck, idcode, needav, level, tol, timsys, room );
   catch
      rethrow(lasterror)
   end

   %
   % Union the coverage window 'cov' returned from the interface to the
   % input coverage window, 'cover_i'.
   %
   cov = cspice_wnunid( cov, cover_i );

end

