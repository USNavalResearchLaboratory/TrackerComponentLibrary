%-Abstract
%
%   CSPICE_SPKPVN returns, for a specified SPK segment and time, the state
%   (position and velocity) of the segment's target body relative to its
%   center of motion.
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
%      handle,
%      descr     respectively, file handle assigned to a SPK file and
%                the descriptor for a segment within the file. Together
%                they determine the ephemeris data from which the state
%                of the body is to be computed.
%
%                [1,1] = size(handle); int32 = class(handle)
%
%                [5,1] = size(descr); double = class(descr)
%
%      et        a scalar double precision time, in seconds past the
%                epoch J2000 TDB.
%
%                [1,1] = size(et); double = class(et)
%
%   the call:
%
%      [ref, state, center] = cspice_spkpvn( handle, descr, et)
%
%   returns:
%
%      ref      the ID code of the reference frame relative to which the
%               state returned by the routine is expressed.
%
%               [1,1] = size(ref); int32 = class(ref)
%
%      state    the array containing the position and velocity, at epoch
%               'et', for the body covered by the specified segment. 'state'
%               has six elements:  the first three contain the body's
%               position; the last three contain the body's velocity. These
%               vectors are expressed into the specified reference frame.
%               Units are always km and km/sec.
%
%               [6,1] = size(state); double = class(state)
%
%      center   the SPK ID code of the center of motion for the state.
%
%               [1,1] = size(center); int32 = class(center)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   In the following code fragment, an SPK file is searched for
%   a segment containing ephemeris data for the Jupiter system
%   barycenter at a particular epoch. Using this segment,
%   states of the Jupiter system barycenter relative to the
%   solar system barycenter are evaluated at a sequence of times.
%
%   This method of state computation minimizes the number of
%   segment searches required to obtain requested data, but
%   it bypasses the SPK subsystem's state chaining mechanism.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%         KPL/MK
%
%         File name: standard.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            de421.bsp                     Planetary ephemeris
%            pck00010.tpc                  Planet orientation and
%                                          radii
%            naif0010.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00010.tpc',
%                                'naif0010.tls'  )
%
%         \begintext
%
%   Example:
%
%      %
%      % Local constants
%      %
%      META   =  'standard.tm';
%      ND     =  2;
%      NI     =  6;
%      TIMFMT =  'YYYY MON DD HR:MN:SC.######::TDB TDB';
%
%      %
%      % Load meta-kernel.
%      %
%      cspice_furnsh( META )
%
%      %
%      % Convert starting time to seconds past J2000 TDB.
%      %
%      timstr = '2012 APR 27 00:00:00.000 TDB';
%
%      et0 = cspice_str2et(timstr);
%
%      %
%      % Find a loaded segment for the Jupiter barycenter
%      % that covers `et0'.
%      %
%      body = 5;
%
%      [handle, descr, segid, found] = cspice_spksfs( body, et0);
%
%
%      if ~found
%         cspice_kclear
%         txt = sprintf( 'No SPK segment found for body %d at time %s', ...
%                         body, timstr );
%         error( txt )
%      end
%
%      %
%      % Unpack the descriptor of the current segment.
%      %
%      [dc, ic] = cspice_dafus( descr, ND, NI );
%
%      frname = cspice_frmnam( ic(3) );
%
%      fprintf( 'Body        = %d\n', ic(1) )
%      fprintf( 'Center      = %d\n', ic(2) )
%      fprintf( 'Frame       = %s\n', frname)
%      fprintf( 'Data type   = %d\n', ic(4) )
%      fprintf( 'Start ET    = %f\n', dc(1) )
%      fprintf( 'Stop ET     = %f\n', dc(2) )
%      fprintf( 'Segment ID  = %s\n\n', segid )
%
%
%      %
%      % Evaluate states at 10-second steps, starting at `et0'
%      % and continuing for 20 seconds.
%      %
%
%      for i=1:3
%
%         et = et0 + ( 10. * (i-1) );
%
%         %
%         % Convert `et' to a string for display.
%         %
%         outstr = cspice_timout( et, TIMFMT );
%
%         %
%         % Attempt to compute a state only if the segment's
%         % coverage interval contains `et'.
%         %
%         if ( et <= dc(2) )
%
%            %
%            % This segment has data at `et'. Evaluate the
%            % state of the target relative to its center
%            % of motion.
%            %
%            [ref_id, state, center] = cspice_spkpvn( handle, descr, et );
%
%            %
%            %  Display the time and state.
%            %
%            fprintf( '\n%s\n', outstr )
%            fprintf( 'Position X (km):   %24.17f\n', state(1) )
%            fprintf( 'Position Y (km):   %24.17f\n', state(2) )
%            fprintf( 'Position Z (km):   %24.17f\n', state(3) )
%            fprintf( 'Velocity X (km):   %24.17f\n', state(4) )
%            fprintf( 'Velocity X (km):   %24.17f\n', state(5) )
%            fprintf( 'Velocity X (km):   %24.17f\n', state(6) )
%
%         else
%
%            cspice_kclear
%            txt = sprintf( 'No data found for body %d at time %s', ...
%                         body, outstr );
%            error( txt )
%
%         end
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in IDL due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Body        = 5
%      Center      = 0
%      Frame       = J2000
%      Data type   = 2
%      Start ET    = -3169195200.000000
%      Stop ET     = 1696852800.000000
%      Segment ID  = DE-0421LE-0421
%
%
%      2012 APR 27 00:00:00.000000 TDB
%      Position X (km):   464528993.98216485977172852
%      Position Y (km):   541513126.15685200691223145
%      Position Z (km):   220785135.62462940812110901
%      Velocity X (km):      -10.38685648307654930
%      Velocity X (km):        7.95324700713742416
%      Velocity X (km):        3.66185835431306517
%
%      2012 APR 27 00:00:10.000000 TDB
%      Position X (km):   464528890.11359262466430664
%      Position Y (km):   541513205.68931341171264648
%      Position Z (km):   220785172.24320945143699646
%      Velocity X (km):      -10.38685796160419272
%      Velocity X (km):        7.95324528430304944
%      Velocity X (km):        3.66185765185608103
%
%      2012 APR 27 00:00:20.000000 TDB
%      Position X (km):   464528786.24500560760498047
%      Position Y (km):   541513285.22175765037536621
%      Position Z (km):   220785208.86178246140480042
%      Velocity X (km):      -10.38685944013147910
%      Velocity X (km):        7.95324356146845002
%      Velocity X (km):        3.66185694939899253
%
%-Particulars
%
%   This routine finds the highest-priority segment, in any loaded
%   SPK file, such that the segment provides data for the specified
%   body and epoch.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine spkpvn_c.
%
%   MICE.REQ
%   SPK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 30-OCT-2012, EDW (JPL)
%
%-Index_Entries
%
%   select spk file and segment
%
%-&

function [ref, state, center] = cspice_spkpvn( handle, descr, et)

   switch nargin
      case 3

         handle = zzmice_int(handle);
         descr  = zzmice_dp(descr);
         et     = zzmice_dp(et);

      otherwise

         error ( ['Usage: [ref, state(6), center] = ' ...
                          'cspice_spkpvn( handle, descr(5), et)'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      spkpvn = mice( 'spkpvn_s', handle, descr, et );

      ref    = reshape( [spkpvn.ref   ], 1, [] );
      state  = reshape( [spkpvn.state ], 6, [] );
      center = reshape( [spkpvn.center], 1, [] );

   catch
      rethrow(lasterror)
   end
