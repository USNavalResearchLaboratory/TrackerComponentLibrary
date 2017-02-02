%-Abstract
%
%   CSPICE_SPKCOV returns the coverage window for a specified ephemeris
%   object in a specified SPK file.
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
%      spk       tthe string name, or cell of strings, of SPICE CK files.
%
%                [1,c] = size(spk), char = class(spk)
%
%                  or
%
%                [1,m] = size(spk), cell = class(spk)
%
%      idcode    the integer scalar NAIF ID code of an object for which
%                ephemeris data are expected to exist in the specified SPK file.
%
%                [1,1] = size(idcode), int32 = class(idcode)
%
%      room      an integer scalar defining the number of intervals for use
%                as a workspace by the routine. This value should equal at least
%                the number of intervals corresponding to 'idcode' in 'spk'.
%
%                [1,1] = size(room), int32 = class(room)
%
%      cover_i   an optional input describing a either an empty window or a
%                window array created from a previous cspice_spkcov call.
%                Inclusion of this window argument results in an output
%                window consisting of a union of the data retrieved from the
%                'spk' kernels and the data in 'cover_i'.
%
%                [2m,1] = size(cover_i), double = class(cover_i)
%
%                   or
%
%                [0,0] = size(cover_i), double = class(cover_i)
%
%   the call:
%
%      cover = cspice_spkcov( spk, idcode, room, cover_i )
%
%         or
%
%      cover = cspice_spkcov( spk, idcode, room )
%
%   returns:
%
%      cover   a double precision 2Nx1 array (window) containing the coverage
%              for 'idcode'', i.e. the set of time intervals for which 'idcode'
%              data exist in the file 'spk'. The array 'cover'
%              contains the pairs of endpoints of these intervals.
%
%               Each window defined as a pair of endpoints such that:
%
%                  window 1 = cover(1:2)
%                  window 2 = cover(3:4)
%                  window 3 = cover(5:6)
%                           ...
%                  window N = cover(2N-1,2N)
%
%              The interval endpoints contained in 'cover' are ephemeris
%              times, expressed as seconds past J2000 TDB.
%
%              'cover' returns an empty set if 'spk' lacks coverage for
%              'idcode'.  If 'cover_i' exists in the argument list, 'cover'
%              returns as a union of the coverage data found in 'spk' and
%              the data in 'cover_i'. 'cover' can overwrite 'cover_i'.
%
%              [2p,1] = size(cover), double = class(cover)
%
%                 or
%
%              [0,1] = size(cover), double = class(cover)
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
%   algorithms using cspice_spkobj will also use cspice_spkcov and
%   vice-versa.
%
%   function spkcov_t(SPK)
%
%      MAXOBJ = 1000;
%      WINSIZ = 1000;
%      LSK    = 'naif0010.tls';
%
%      %
%      % Note, neither cspice_spkcov or cspice_spkobj requires this
%      % kernel to function. We need the data for output time
%      % conversion.
%      %
%      cspice_furnsh( LSK )
%
%      %
%      % Find the set of objects in the SPK file.
%      %
%      ids = cspice_spkobj( SPK, MAXIV );
%
%      %
%      % We want to display the coverage for each object. Loop over
%      % the contents of the ID code set, find the coverage for
%      % each item in the set, and display the coverage.
%      %
%      for i=1:numel(ids)
%
%         %
%         % Extract the coverage data for object 'ids(i)'.
%         %
%         cover     = cspice_spkcov( SPK, ids(i), WINSIZ );
%         [row,col] = size(cover);
%
%         %
%         % Display a simple banner.
%         %
%         fprintf( '========================================\n')
%         fprintf( 'Coverage for object %d\n', ids(i) )
%
%         %
%         %  'cover' has dimension 2Nx1, where 'row' has the value 2N with
%         %  each window defined as a pair of endpoints such that:
%         %
%         %  window 1 = cover(1:2)
%         %  window 2 = cover(3:4)
%         %  window 3 = cover(5:6)
%         %        ...
%         %  window N = cover(2N-1,2N)
%         %
%         % Loop from 1 to 'row' with step size 2.
%         %
%         for j=1:2:row
%
%            %
%            % Convert the endpoints to TDB calendar format time strings
%            % and display them. Pass the endpoints in an array,
%            % so cspice_timout returns an array of time strings.
%            %
%            % Recall a vectorized input has dimension 1xM so transpose
%            % the 'cover' slice.
%            %
%            timstr = cspice_timout( cover(j:j+1)', ...
%                                'YYYY MON DD HR:MN:SC.### (TDB) ::TDB' );
%            fprintf('Interval: %d\n'  , (j+1)/2 )
%            fprintf('   Start: %s\n'  , timstr(1,:) )
%            fprintf('    Stop: %s\n\n', timstr(2,:) )
%
%         end
%
%      end
%
%      %
%      % Empty the kernel pool.
%      %
%      cspice_kclear
%
%   Example (1):
%
%      Assign an SPK kernel list as:
%
%      >> SPK = { '/kernels/gen/spk/de405_2000-2050.bsp', ...
%                 '/kernels/gen/spk/jup100.bsp' };
%
%      >> spkobj_t(SPK)
%
%   MATLAB outputs:
%
%      ========================================
%      Coverage for object 1
%      Interval: 1
%         Start: 2000 JAN 01 00:01:04.183 (TDB)
%          Stop: 2050 JAN 01 00:01:04.183 (TDB)
%
%      ========================================
%      Coverage for object 2
%      Interval: 1
%         Start: 2000 JAN 01 00:01:04.183 (TDB)
%          Stop: 2050 JAN 01 00:01:04.183 (TDB)
%
%      ========================================
%      Coverage for object 3
%      Interval: 1
%         Start: 1973 NOV 01 00:00:00.000 (TDB)
%          Stop: 2050 JAN 01 00:01:04.183 (TDB)
%
%      ========================================
%      Coverage for object 4
%      Interval: 1
%         Start: 2000 JAN 01 00:01:04.183 (TDB)
%          Stop: 2050 JAN 01 00:01:04.183 (TDB)
%
%      ========================================
%      Coverage for object 5
%      Interval: 1
%         Start: 1973 NOV 01 00:00:00.000 (TDB)
%          Stop: 2050 JAN 01 00:01:04.183 (TDB)
%
%             ... continued ...
%
%      ========================================
%      Coverage for object 501
%      Interval: 1
%         Start: 1973 NOV 01 00:00:00.000 (TDB)
%          Stop: 2023 NOV 03 00:00:00.000 (TDB)
%
%      ========================================
%      Coverage for object 502
%      Interval: 1
%         Start: 1973 NOV 01 00:00:00.000 (TDB)
%          Stop: 2023 NOV 03 00:00:00.000 (TDB)
%
%      ========================================
%      Coverage for object 503
%      Interval: 1
%         Start: 1973 NOV 01 00:00:00.000 (TDB)
%          Stop: 2023 NOV 03 00:00:00.000 (TDB)
%
%      ========================================
%      Coverage for object 504
%      Interval: 1
%         Start: 1973 NOV 01 00:00:00.000 (TDB)
%          Stop: 2023 NOV 03 00:00:00.000 (TDB)
%
%      ========================================
%      Coverage for object 599
%      Interval: 1
%         Start: 1973 NOV 01 00:00:00.000 (TDB)
%          Stop: 2023 NOV 03 00:00:00.000 (TDB)
%
%   Example (2):
%
%      Assign the SPK kernel list as:
%
%      >> SPK = { '/kernels/Hubble/1990-01-01_1996-01-01.bsp', ...
%                 '/kernels/Hubble/2002-01-01_2006-12-31.bsp' };
%
%      >> spkobj_t(SPK)
%
%   MATLAB outputs:
%
%      ========================================
%      Coverage for object -48
%      Interval: 1
%         Start: 1990 APR 25 02:01:40.071 (TDB)
%          Stop: 1995 DEC 31 17:42:03.755 (TDB)
%
%      Interval: 2
%         Start: 2002 JAN 01 02:02:00.084 (TDB)
%          Stop: 2006 DEC 29 09:36:39.732 (TDB)
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine spkcov_c.
%
%   MICE.REQ
%   SPK.REQ
%   CELLS.REQ
%   DAF.REQ
%   TIME.REQ
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
%      Corrected several typos.
%
%   -Mice Version 1.0.0, 18-JUN-2007, EDW (JPL)
%
%-Index_Entries
%
%   get coverage window for spk object
%
%-&

function [cover] = cspice_spkcov( spk, idcode, room, cover_i )

   switch nargin
      case 3

         spk    = zzmice_str(spk);
         idcode = zzmice_int(idcode);
         room   = zzmice_int(room, [1, int32(inf)/2] );

      case 4

         spk    = zzmice_str(spk);
         idcode = zzmice_int(idcode);
         room   = zzmice_int(room, [1, int32(inf)/2] );
         cover_i= zzmice_win(cover_i);

      otherwise

         error ( [ 'Usage: [cover] = cspice_spkcov( _`spk`_, ' ...
                                     'idcode, room, [cover_i])' ] )

   end

%
% The call passed either three or four arguments. Branch accordingly.
%
if ( nargin == 3 )

   %
   % Call the MEX library.
   %
   try
      [cover] = mice('spkcov_c', spk, idcode, room );
   catch
      rethrow(lasterror)
   end

else

   %
   % Call the MEX library.
   %
   try
      cover = mice('spkcov_c', spk, idcode, room );
   catch
      rethrow(lasterror)
   end

   cover = cspice_wnunid( cover, cover_i );

end



