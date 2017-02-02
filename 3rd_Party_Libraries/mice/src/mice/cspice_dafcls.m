%-Abstract
%
%   CSPICE_DAFCLS closes the DAF referred to by a handle.
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
%      handle   file handle referring to a DAF.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      cspice_dafcls( handle )
%
%   closes the DAF, i.e. removes the file from read or write
%   access by the SPICE DAF subsystem.
%
%   Use this routine to close files opened by cspice_dafopr and
%   cspice_dafopw.
%
%   returns:
%
%      None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Use a simple function to output the double precision and integer
%   values stored in an SPK's segments descriptors. This function opens
%   a DAF for read, performs a forwards search for the DAF arrays,
%   prints the segment description for each array found, then closes the DAF.
%
%   function daf_t( kernel)
%
%      %
%      % Open a DAF for read. Return a 'handle' referring to the file.
%      %
%      handle = cspice_dafopr( kernel );
%
%      %
%      % Define the summary parameters appropriate
%      % for an SPK file.
%      %
%      ND = 2;
%      NI = 6;
%
%      %
%      % Begin a forward search on the file.
%      %
%      cspice_dafbfs( handle );
%
%      %
%      % Search until a DAF array is found.
%      %
%      found = cspice_daffna;
%
%      %
%      % Loop while the search finds subsequent DAF arrays.
%      %
%      while found
%
%         [dc, ic ] = cspice_dafgs( ND, NI );
%
%         fprintf( 'Doubles:  ' )
%         fprintf( '%f   ', dc )
%         fprintf( '\n' )
%
%         fprintf( 'Integers: ' )
%         fprintf( '%d   ', ic )
%         fprintf( '\n\n' )
%
%
%         %
%         % Check for another segment.
%         %
%         found = cspice_daffna;
%
%      end
%
%      %
%      % Safely close the DAF.
%      %
%      cspice_dafcls( handle )
%
%   Matlab outputs:
%
%      >> daf_t( 'de421.bsp' )
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 1   0   1   2   641   310404
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 2   0   1   2   310405   423048
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 3   0   1   2   423049   567372
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 4   0   1   2   567373   628976
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 5   0   1   2   628977   674740
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 6   0   1   2   674741   715224
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 7   0   1   2   715225   750428
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 8   0   1   2   750429   785632
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 9   0   1   2   785633   820836
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 10   0   1   2   820837   944040
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 301   3   1   2   944041   1521324
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 399   3   1   2   1521325   2098608
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 199   1   1   2   2098609   2098620
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 299   2   1   2   2098621   2098632
%
%      Doubles:  -3169195200.000000   1696852800.000000
%      Integers: 499   4   1   2   2098633   2098644
%
%   Note, the specific contents of 'ic' and 'dc' depend on the
%   type of DAF.
%
%   Note, the final entries in the integer array contain the segment
%   start/end indexes. The output indicates the search proceeded
%   from the start of the file (low value index) towards the end
%   (high value index).
%
%-Particulars
%
%   Because the DAF subsystem must keep track of what files are open at
%   any given time, it is important that DAF files be closed only with
%   cspice_dafcls, to prevent the remaining DAF routines from failing,
%   sometimes mysteriously.
%
%   Note that when a file is opened more than once for read access,
%   cspice_dafopr returns the same handle each time it is re-opened.
%   Each time the file is closed, cspice_dafcls checks to see if any other
%   claims on the file are still active before physically closing
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dafcls_c.
%
%   MICE.REQ
%   DAF.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 10-JUL-2012, EDW (JPL)
%
%-Index_Entries
%
%   close DAF
%
%-&

function cspice_dafcls( handle )

   switch nargin
      case 1

         handle  = zzmice_int(handle);

      otherwise

         error ( 'Usage: cspice_dafcls(handle)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'dafcls_c', handle );
   catch
      rethrow(lasterror)
   end




