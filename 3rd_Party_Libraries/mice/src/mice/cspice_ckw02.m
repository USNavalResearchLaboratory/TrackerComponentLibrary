%-Abstract
%
%   CSPICE_CKW02 adds a type 2 segment to a C-kernel.
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
%      handle   scalar integer value of the file handle for
%               an open CK file returned from cspice_ckopn
%
%      begtim   double precision scalar encoded SCLK segment
%               begin time
%
%      endtim   double precision scalar encoded SCLK segment
%               end time
%
%      inst     the scalar integer NAIF instrument ID code
%
%      ref      scalar string identifying the reference frame for the
%               segment
%
%      segid    a scalar string to identify the segment
%
%      start    a double precision Nx1 array containing
%               encoded SCLK interval start times
%
%      stop     a double precision Nx1 array containing
%               the encoded SCLK interval stop times
%
%      quats    a double precision 4xN matrix of SPICE style quaternions
%               representing instrument pointing
%
%      avvs     a double precision 3xN  matrix of angular
%               velocity vectors in units of radians per second
%
%      rates    a double precision Nx1 array containing the
%               number of seconds per tick for each interval
%
%   the call:
%
%      cspice_ckw02( handle,  ...
%                    begtime, ...
%                    endtime, ...
%                    inst,    ...
%                    ref,     ...
%                    segid,   ...
%                    start,   ...
%                    stop,    ...
%                    quats,   ...
%                    avvs,    ...
%                    rates )
%
%   returns:
%
%      Adds a type 2 segment to a CK.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Create a CK type 2 segment; fill with data for a simple time
%      % dependent rotation and angular velocity.
%      %
%
%      %
%      % Define needed parameters
%      %
%      CK2        = 'type2.bc';
%      INST       = -77702;
%      MAXREC     = 21;
%      SECPERTICK = 0.001;
%      IFNAME     = 'Test CK type 2 segment created by cspice_ckw02';
%      SEGID      = 'Test type 2 CK segment';
%
%      %
%      % 'NCOMCH' is the number of characters to reserve for the kernel's
%      % comment area. This example doesn't write comments, so set to
%      % zero.
%      %
%      NCOMCH = 0;
%
%      %
%      % The base reference from for the rotation data.
%      %
%      REF = 'J2000';
%
%      %
%      % Time spacing in encoded ticks.
%      %
%      SPACING_TICKS = 10.;
%
%      %
%      % Time spacing in seconds
%      %
%      SPACING_SECS = SPACING_TICKS * SECPERTICK;
%
%      %
%      % Declare an angular rate in radians per sec.
%      %
%      RATE = 1.d-2;
%
%      %
%      % Create a 4xMAXREC matrix for quaternions, and a
%      % 3xMAXREC for expavs.
%      %
%      quats = zeros( [4, MAXREC] );
%      av    = zeros( [3, MAXREC] );
%
%      %
%      % Create a 3x3 double precision identity matrix.
%      %
%      work_mat = eye( 3 );
%
%      %
%      % Convert the 'work_mat' to quaternions.
%      %
%      work_quat = cspice_m2q( work_mat);
%
%      %
%      % Copy the work quaternions to the first row of
%      % 'quats'.
%      %
%      quats(:,1) = work_quat;
%
%      %
%      % Create an angular velocity vector. Copy to the third (Z) row
%      % of 'av'. This vector is in the 'REF' reference frame and
%      % indicates a constant rotation about the Z axis.
%      %
%      av(3,:) = RATE;
%
%      %
%      % Create arrays of interval start and stop times.  The interval
%      % associated with each quaternion will start at the epoch of
%      % the quaternion and will extend 0.8 * SPACING_TICKS forward in time,
%      % leaving small gaps between the intervals.
%      %
%      % Fill in the clock rates array with a constant 'SECPERTICK' for
%      % all values.
%      %
%      rates  = zeros( [MAXREC,1] ) + SECPERTICK;
%
%      %
%      % Create an array of encoded tick values in increments of
%      % 'SPACING_TICKS' with an initial value of 1000 ticks...
%      %
%      sclkdp = [0:MAXREC-1]' * SPACING_TICKS + 1000;
%
%      starts = sclkdp;
%      stops  = sclkdp + ( 0.8 * SPACING_TICKS );
%
%      %
%      % Fill the rest of the av and quats matrices
%      % with simple data.
%      %
%      for i = 2:MAXREC
%
%         %
%         % Create the transformation matrix for a rotation of 'theta'
%         % about the Z axis. Calculate 'theta' from the constant
%         % angular rate 'RATE' at increments of 'SPACING_SECS'.
%         %
%         theta    = (i-1) * RATE * SPACING_SECS;
%         work_mat = cspice_rotmat( work_mat, theta, 3 );
%
%         %
%         % Convert the 'work_mat' matrix to SPICE type quaternions.
%         %
%         work_quat = cspice_m2q( work_mat );
%
%         %
%         % Store the quaternions in the 'quat' matrix.
%         %
%         quats(:,i) = work_quat;
%
%      end
%
%      %
%      % Set the segment boundaries equal to the first and last
%      % time in the segment.
%      %
%      begtime = starts(1);
%      endtime = stops(MAXREC);
%
%      %
%      % All information ready to write. Write to a CK type 2 segment
%      % to the file indicated by 'handle'.
%      %
%      try
%         handle = cspice_ckopn( CK2, IFNAME, NCOMCH );
%         cspice_ckw02(  handle,  ...
%                        begtime, ...
%                        endtime, ...
%                        INST,    ...
%                        REF,     ...
%                        SEGID,   ...
%                        starts,  ...
%                        stops,   ...
%                        quats,   ...
%                        av,      ...
%                        rates )
%      catch
%         error( [ 'Failure: ' lasterr] )
%      end
%
%      %
%      % SAFELY close the file.
%      %
%      cspice_ckcls( handle )
%
%   MATLAB outputs:
%
%      The example code creates a CK with one type 2 segment.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ckw02_c.
%
%   MICE.REQ
%   CK.REQ
%   DAF.REQ
%   SCLK.REQ
%
%-Version
%
%   -Mice Version 1.0.2, 11-JUN-2013, EDW (JPL)
%
%       I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.1, 30-DEC-2008, EDW (JPL)
%
%      Corrected misspellings.
%
%   -Mice Version 1.0.0, 04-JAN-2008, EDW (JPL)
%
%-Index_Entries
%
%   write ck type_2 pointing data segment
%
%-&

function cspice_ckw02( handle,  ...
                       begtime, ...
                       endtime, ...
                       inst,    ...
                       ref,     ...
                       segid,   ...
                       start,   ...
                       stop,    ...
                       quats,   ...
                       avvs,    ...
                       rates )

   switch nargin
      case 11

         handle  = zzmice_int(handle);
         begtime = zzmice_dp(begtime);
         endtime = zzmice_dp(endtime);
         inst    = zzmice_int(inst);
         ref     = zzmice_str(ref);
         segid   = zzmice_str(segid);
         start   = zzmice_dp(start);
         stop    = zzmice_dp(stop);
         quats   = zzmice_dp(quats);
         avvs    = zzmice_dp(avvs);
         rates   = zzmice_dp(rates);

      otherwise

         error ( [ 'Usage: '                    ...
                   'cspice_ckw02( handle, '     ...
                                 'begtime, '    ...
                                 'endtime, '    ...
                                 'inst, '       ...
                                 'ref, '        ...
                                 'segid, '      ...
                                 'start(N), '   ...
                                 'stop(N), '    ...
                                 'quats(4,N), ' ...
                                 'avvs(3,N)) '  ...
                                 'rates(N)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'ckw02_c', handle,  ...
                       begtime, ...
                       endtime, ...
                       inst,    ...
                       ref,     ...
                       segid,   ...
                       start,   ...
                       stop,    ...
                       quats,   ...
                       avvs,    ...
                       rates )
   catch
      rethrow(lasterror)
   end


