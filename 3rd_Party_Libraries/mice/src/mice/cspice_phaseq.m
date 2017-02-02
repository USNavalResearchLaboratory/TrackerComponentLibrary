%-Abstract
%
%   CSPICE_PHASEQ computes the apparent phase angle for a target, observer,
%   illuminator set of ephemeris objects.
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
%      et       the epochs, specified in ephemeris seconds past J2000, at which
%               to compute the phase angle.
%
%               [1,n] = size(et), double = class(et)
%
%      target   the string naming of the target body.
%
%               Optionally, you may supply the integer NAIF ID code
%               for the body as a string. For example both 'MOON' and
%               '301' are legitimate strings that designate the Moon.
%
%               Case and leading or trailing blanks are not significant
%               in the string 'target'.
%
%               [1,c1] = size(target), char = class(target)
%
%      illmn    the string naming the illuminating body.
%
%               Optionally, you may supply the integer NAIF ID code
%               for the body as a string. For example both 'MOON' and
%               '301' are legitimate strings that designate the Moon.
%
%               Case and leading or trailing blanks are not significant
%               in the string 'illmn'.
%
%               In most cases, 'illmn' is the sun.
%
%               [1,c2] = size(target), char = class(target)
%
%      obsrvr   the string naming the observing body, typically a
%               spacecraft, the earth, or a surface point on the earth.
%
%               Optionally, you may supply the integer NAIF ID code
%               for the body as a string. For example both 'MOON' and
%               '301' are legitimate strings that designate the Moon.
%
%               Case and leading or trailing blanks are not significant
%               in the string 'obsrvr'.
%
%               [1,c3] = size(obsrvr), char = class(obsrvr)
%
%      abcorr   the string naming the aberration corrections to apply
%               to the state evaluations to account for one-way light time and
%               stellar aberration.
%
%               This routine accepts only reception mode aberration
%               corrections. See the header of cspice_spkezr for a detailed
%               description of the aberration correction options.
%               For convenience, the appropriate aberration options are
%               listed below:
%
%                  'NONE'     Apply no correction. Returns the "true"
%                             geometric state.
%
%                  'LT'       "Reception" case:  correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'LT+S'     "Reception" case:  correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'CN'       "Reception" case:  converged
%                             Newtonian light time correction.
%
%                  'CN+S'     "Reception" case:  converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%               Case and leading or trailing blanks are not significant
%               in the string 'abcorr'.
%
%               [1,c4] = size(abcorr), char = class(abcorr)
%
%   the call:
%
%      phase = cspice_phaseq( et, target, illum, obsrvr, abcorr )
%
%   returns:
%
%      phase   the optionally light-time corrected phase angle between
%              'target' and 'illmn' as observed  from 'obsrvr'.
%              Units are radians.  The range of 'phase' is [0, pi].
%
%              'phase' return with the same vectorization measure (N) as 'et'.
%
%              [1,n] = size(phase), double = class(phase)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
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
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0009.tls'
%                                'de421.bsp'
%                                'pck00009.tpc' )
%
%         \begintext
%
%   Example:
%
%      Determine the time intervals from December 1, 2006 UTC to
%      January 31, 2007 UTC for which the sun-moon-earth configuration
%      phase angle satisfies the relation conditions with respect to a
%      reference value of .57598845 radians (the phase angle at
%      January 1, 2007 00:00:00.000 UTC, 33.001707 degrees). Also
%      determine the time intervals corresponding to the local maximum and
%      minimum phase angles, and the absolute maximum and minimum phase
%      angles during the search interval. The configuration defines the
%      sun as the illuminator, the moon as the target, and the earth as
%      the observer.
%
%      MAXWIN  =  5000;
%      TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###';
%
%      relate = { '=', '<', '>', ...
%                 'LOCMIN', 'ABSMIN', 'LOCMAX', 'ABSMAX' };
%
%      %
%      % Define the location for the phase angle calculation as the
%      % geometric center of the target.
%      %
%      pos = [ 0, 0, 0 ]';
%
%      %
%      % Load kernels.
%      %
%      cspice_furnsh( 'standard.tm' );
%
%      %
%      % Store the time bounds of our search interval in
%      % the cnfine confinement window.
%      %
%      et = cspice_str2et( { '2006 DEC 01', '2007 JAN 31'} );
%
%      %
%      % Search using a step size of 1 day (in units of seconds).
%      % The reference value is 0.57598845 radians. We're not using the
%      % adjustment feature, so we set 'adjust' to zero.
%      %
%      target  = 'MOON';
%      illum   = 'SUN';
%      abcorr  = 'LT+S';
%      obsrvr  = 'EARTH';
%      refval  = 0.57598845;
%      adjust  = 0.;
%      step    = cspice_spd;
%      nintvls = MAXWIN;
%      cnfine  = cspice_wninsd( et(1), et(2) );
%
%      for j=1:numel( relate )
%
%         fprintf( 'Relation condition: %s\n',  char( relate(j) ) )
%
%         %
%         % Perform the search. The SPICE window 'result' contains
%         % the set of times when the condition is met.
%         %
%         result = cspice_gfpa( target,    illum,  abcorr, obsrvr, ...
%                               relate(j), refval, adjust, step,  ...
%                               nintvls,   cnfine );
%
%         %
%         % Display the results.
%         %
%         count = cspice_wncard(result);
%
%         if ( isequal( count, 0 ) )
%
%               fprintf( 'Result window is empty.\n\n' );
%
%         else
%
%            for i=1:count
%
%               %
%               % Fetch the endpoints of the Ith interval
%               % of the result window.
%               %
%               [left, right] = cspice_wnfetd( result, i );
%
%               phase = cspice_phaseq( [left, right], target, illum, ...
%                                      obsrvr, abcorr );
%
%               output = cspice_timout( [left,right], TIMFMT );
%
%               fprintf( 'Start time = %s %16.9f\n', output(1,:), phase(1) )
%               fprintf( 'Stop time  = %s %16.9f\n', output(2,:), phase(2) )
%
%            end
%
%            disp( ' ')
%
%         end
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Relation condition: =
%      Start time = 2006-DEC-02 13:31:34.414      0.575988450
%      Stop time  = 2006-DEC-02 13:31:34.414      0.575988450
%      Start time = 2006-DEC-07 14:07:55.470      0.575988450
%      Stop time  = 2006-DEC-07 14:07:55.470      0.575988450
%      Start time = 2006-DEC-31 23:59:59.997      0.575988450
%      Stop time  = 2006-DEC-31 23:59:59.997      0.575988450
%      Start time = 2007-JAN-06 08:16:25.512      0.575988450
%      Stop time  = 2007-JAN-06 08:16:25.512      0.575988450
%      Start time = 2007-JAN-30 11:41:32.557      0.575988450
%      Stop time  = 2007-JAN-30 11:41:32.557      0.575988450
%
%      Relation condition: <
%      Start time = 2006-DEC-02 13:31:34.414      0.575988450
%      Stop time  = 2006-DEC-07 14:07:55.470      0.575988450
%      Start time = 2006-DEC-31 23:59:59.997      0.575988450
%      Stop time  = 2007-JAN-06 08:16:25.512      0.575988450
%      Start time = 2007-JAN-30 11:41:32.557      0.575988450
%      Stop time  = 2007-JAN-31 00:00:00.000      0.468279091
%
%      Relation condition: >
%      Start time = 2006-DEC-01 00:00:00.000      0.940714974
%      Stop time  = 2006-DEC-02 13:31:34.414      0.575988450
%      Start time = 2006-DEC-07 14:07:55.470      0.575988450
%      Stop time  = 2006-DEC-31 23:59:59.997      0.575988450
%      Start time = 2007-JAN-06 08:16:25.512      0.575988450
%      Stop time  = 2007-JAN-30 11:41:32.557      0.575988450
%
%      Relation condition: LOCMIN
%      Start time = 2006-DEC-05 00:16:50.317      0.086121423
%      Stop time  = 2006-DEC-05 00:16:50.317      0.086121423
%      Start time = 2007-JAN-03 14:18:31.977      0.079899769
%      Stop time  = 2007-JAN-03 14:18:31.977      0.079899769
%
%      Relation condition: ABSMIN
%      Start time = 2007-JAN-03 14:18:31.977      0.079899769
%      Stop time  = 2007-JAN-03 14:18:31.977      0.079899769
%
%      Relation condition: LOCMAX
%      Start time = 2006-DEC-20 14:09:10.392      3.055062862
%      Stop time  = 2006-DEC-20 14:09:10.392      3.055062862
%      Start time = 2007-JAN-19 04:27:54.600      3.074603891
%      Stop time  = 2007-JAN-19 04:27:54.600      3.074603891
%
%      Relation condition: ABSMAX
%      Start time = 2007-JAN-19 04:27:54.600      3.074603891
%      Stop time  = 2007-JAN-19 04:27:54.600      3.074603891
%
%-Particulars
%
%   This routine returns the phase angle using the location of the
%   bodies (if point objects) or the centers of the bodies (if finite
%   bodies).
%
%
%                       ILLUM      OBS
%       ILLUM as seen      ^       /
%       from TARG at       |      /
%       ET - LT.           |     /
%                         >|..../< phase angle
%                          |   /
%                        . |  /
%                      .   | /
%                     .    |v     TARG as seen from OBS
%               SEP   .   TARG    at ET
%                      .  /
%                        /
%                       v
%
%        PI = SEP + PHASE
%
%        so
%
%        PHASE = PI - SEP
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine phaseq_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 13-MAR-2012, EDW (JPL)
%
%-Index_Entries
%
%   compute the phase of two objects wrt an illumination source
%
%-&

function [phase] = cspice_phaseq( et, target, illmn, obsrvr, abcorr )

   switch nargin
      case 5

         et     = zzmice_dp(et);
         target = zzmice_str(target);
         illmn = zzmice_str(illmn);
         obsrvr = zzmice_str(obsrvr);
         abcorr = zzmice_str(abcorr);

      otherwise

         error ( ['Usage: [_phase_] = cspice_phaseq( _et_, '       ...
                                  '`target`, `illmn`, `obsrvr`, `abcorr` )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [phase] = mice('phaseq_c', et, target, illmn, obsrvr, abcorr );
   catch
      rethrow(lasterror)
   end




