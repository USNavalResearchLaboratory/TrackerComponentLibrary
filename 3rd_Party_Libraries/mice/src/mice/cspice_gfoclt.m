%-Abstract
%
%   CSPICE_GFOCLT determines time intervals when an observer sees one target
%   body occulted by, or in transit across, another.
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
%      Parameters-
%
%      All parameters described here are declared in the header file
%      SpiceGF.h. See that file for parameter values.
%
%      SPICE_GF_CNVTOL
%
%               is the convergence tolerance used for finding endpoints of
%               the intervals comprising the result window.
%               SPICE_GF_CNVTOL is used to determine when binary searches
%               for roots should terminate: when a root is bracketed
%               within an interval of length SPICE_GF_CNVTOL, the root is
%               considered to have been found.
%
%               The accuracy, as opposed to precision, of roots found
%               by this routine depends on the accuracy of the input
%               data. In most cases, the accuracy of solutions will be
%               inferior to their precision.
%
%      Arguments-
%
%      occtyp   the string naming the type of occultation to find.
%               Note that transits are considered to be a type of
%               occultation.
%
%               Supported values and corresponding definitions are:
%
%                  'FULL'               denotes the full occultation
%                                       of the body designated by
%                                       'back' by the body designated
%                                       by 'front', as seen from
%                                       the location of the observer.
%                                       In other words, the occulted
%                                       body is completely invisible
%                                       as seen from the observer's
%                                       location.
%
%                  'ANNULAR'            denotes an annular
%                                       occultation: the body
%                                       designated by 'front' blocks
%                                       part of, but not the limb of,
%                                       the body designated by 'back',
%                                       as seen from the location of
%                                       the observer.
%
%                  'PARTIAL'            denotes a partial,
%                                       non-annular occultation: the
%                                       body designated by 'front'
%                                       blocks part, but not all, of
%                                       the limb of the body
%                                       designated by 'back', as seen
%                                       from the location of the
%                                       observer.
%
%                  'ANY'                denotes any of the above three
%                                       types of occultations:
%                                       'PARTIAL', 'ANNULAR', or
%                                       'FULL'.
%
%                                       'ANY' should be used to search
%                                       for times when the body
%                                       designated by 'front' blocks
%                                       any part of the body designated
%                                       by 'back'.
%
%                                       The option 'ANY' must be used
%                                       if either the front or back
%                                       target body is modeled as
%                                       a point.
%
%               Case and leading or trailing blanks are not significant in
%               the string 'occtyp'.
%
%               [1,c1] = size(occtyp), char = class(occtyp)
%
%      front    the string naming the target body that occults---that
%               is, passes in front of---the other. Optionally, you may
%               supply the integer NAIF ID code for the body as a
%               string. For example both 'MOON' and '301' are
%               legitimate strings that designate the Moon.
%
%               The 'front' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%               [1,c2] = size(front), char = class(front)
%
%      fshape   the string naming the geometric model used
%               to represent the shape of the front target body. The
%               supported options are:
%
%                 'ELLIPSOID'     Use a triaxial ellipsoid model
%                                 with radius values provided from the
%                                 kernel pool. A kernel variable
%                                 having a name of the form
%
%                                    'BODYnnn_RADII'
%
%                                 where nnn represents the NAIF
%                                 integer code associated with the
%                                 body, must be present in the kernel
%                                 pool. This variable must be
%                                 associated with three numeric
%                                 values giving the lengths of the
%                                 ellipsoid's X, Y, and Z semi-axes.
%
%                 'POINT'         Treat the body as a single point.
%                                 When a point target is specified,
%                                 the occultation type must be
%                                 set to 'ANY'.
%
%               At least one of the target bodies 'front' and 'back' must
%               be modeled as an ellipsoid.
%
%               The 'fshape' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%               [1,c3] = size(fshape), char = class(fshape)
%
%      fframe   the string naming the body-fixed, body-centered reference
%               frame associated with the front target body. Examples
%               of such names are 'IAU_SATURN' (for Saturn) and
%               'ITRF93' (for the Earth).
%
%               If the front target body is modeled as a point, 'fframe'
%               should be left empty or blank.
%
%               The 'fframe' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%               [1,c4] = size(fframe), char = class(fframe)
%
%      back     the string naming the target body that is occulted
%               by---that is, passes in back of---the other.
%               Optionally, you may supply the integer NAIF ID code
%               for the body as a string. For example both 'MOON' and
%               '301' are legitimate strings that designate the Moon.
%
%               The 'back' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%               [1,c5] = size(back), char = class(back)
%
%      bshape   the string naming the shape specification for the body
%               designated by 'back'. The supported options are those for
%                'fshape'. See the description of 'fshape' above for
%                details.
%
%               [1,c6] = size(bshape), char = class(bshape)
%
%      bframe   the string naming the body-fixed, body-centered
%               reference frame associated with the ''back'' target body.
%               Examples of such names are 'IAU_SATURN' (for Saturn)
%               and 'ITRF93' (for the Earth).
%
%               If the back target body is modeled as a point, 'bframe'
%               should be left empty or blank.
%
%               The 'bframe' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%               [1,c7] = size(bframe), char = class(bframe)
%
%      abcorr   the string indicating the aberration corrections to to apply
%               to the state of the target body to account for one-way
%               light time.  Stellar aberration corrections are
%               ignored if specified, since these corrections don't
%               improve the accuracy of the occultation determination.
%
%               This routine accepts the same aberration corrections as does
%               the CSPICE routine spkezr_c. See the header of spkezr_c for a
%               detailed description of the aberration correction options.
%               For convenience, the options are listed below:
%
%                  'NONE'     Apply no correction.
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
%                  'XLT'      "Transmission" case:  correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'XLT+S'    "Transmission" case:  correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'XCN'      "Transmission" case:  converged
%                             Newtonian light time correction.
%
%                  'XCN+S'    "Transmission" case:  converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%               The 'abcorr' string lacks sensitivity to case, and to embedded,
%               leading and trailing blanks.
%
%               [1,c8] = size(abcorr), char = class(abcorr)
%
%      obsrvr   the scalar naming the observing body. Optionally, you
%               may supply the ID code of the object as an integer string.
%               For example, both 'EARTH' and '399' are legitimate
%               strings to supply to indicate the observer is Earth.
%
%               Case and leading or trailing blanks are not significant in
%               the string 'obsrvr'.
%
%               [1,c9] = size(obsrvr), char = class(obsrvr)
%
%      step     the step size to use in the search. 'step' must be shorter
%               than any interval, within the confinement window, over which
%               the specified occultation condition is met. In other words,
%               'step' must be shorter than the shortest occultation event
%               the user wishes to detect; 'step' must also be shorter than
%               the shortest time interval between two occultation events that
%               occur within the confinement window (see below). However,
%               'step' must not be *too* short, or the search will take
%               an unreasonable amount of time.
%
%               The choice of 'step' affects the completeness but not
%               the precision of solutions found by this routine; the
%               precision is controlled by the convergence tolerance.
%               See the discussion of the parameter SPICE_GF_CNVTOL for
%               details.
%
%               'step' has units of TDB seconds.
%
%               [1,1] = size(step), double = class(step)
%
%      cnfine   the SPICE window that confines the time
%               period over which the specified search is conducted.
%               'cnfine' may consist of a single interval or a collection
%               of intervals.
%
%               In some cases the confinement window can be used to
%               greatly reduce the time period that must be searched
%               for the desired solution. See the Particulars section
%               below for further discussion.
%
%               [2m,1] = size(cnfine), double = class(cnfine)
%
%      room     the maximum number of intervals to return in 'result'.
%               Note: this value should equal at least the number of expected
%               intervals. Recall two double precision values define
%               an interval.
%
%               [1,1] = size(room), int32 = class(room)
%
%   the call:
%
%      result = cspice_gfoclt( occtyp, front,  fshape, fframe, ...
%                              back,   bshape, bframe, abcorr, ...
%                              obsrvr, step,   cnfine, room)
%
%   returns:
%
%      result   the SPICE window of intervals, contained within the
%               confinement window 'cnfine', on which the specified
%               constraint is satisfied.
%
%               If no times within the confinement window satisfy the
%               constraint, 'result' will return with cardinality zero.
%
%               [2n,1] = size(result), double = class(result)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      Find occultations of the Sun by the Moon (that is, solar
%      eclipses)as seen from the center of the Earth over the month
%      December, 2001.
%
%      Use light time corrections to model apparent positions of Sun
%      and Moon. Stellar aberration corrections are not specified
%      because they don't affect occultation computations.
%
%      We select a step size of 3 minutes, which means we
%      ignore occultation events lasting less than 3 minutes,
%      if any exist.
%
%
%      MAXWIN  =  1000;
%      TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
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
%      et = cspice_str2et( { '2001 DEC 01 00:00:00 TDB', ...
%                            '2002 JAN 01 00:00:00 TDB'} );
%
%      cnfine = cspice_wninsd( et(1), et(2) );
%
%      %
%      % Select a 3-minute step. We'll ignore any occultations
%      % lasting less than 3 minutes.
%      %
%      step    = 180.;
%
%      occtyp  = 'any';
%      front   = 'moon';
%      fshape  = 'ellipsoid';
%      fframe  = 'iau_moon';
%      back    = 'sun';
%      bshape  = 'ellipsoid';
%      bframe  = 'iau_sun';
%      obsrvr  = 'earth';
%      abcorr  = 'lt';
%
%      result = cspice_gfoclt( occtyp, front, fshape, fframe, ...
%                              back, bshape, bframe,          ...
%                              abcorr, obsrvr, step, cnfine,  ...
%                              MAXWIN);
%
%      %
%      % List the beginning and ending times in each interval
%      % if result contains data.
%      %
%      for i=1:numel(result)/2
%
%         [left, right] = cspice_wnfetd( result, i );
%
%         output = cspice_timout( [left,right], TIMFMT );
%
%         if( isequal( left, right) )
%
%            disp( ['Event time: ' output(1,:)] )
%
%         else
%
%            disp( ['From : ' output(1,:)] )
%            disp( ['To   : ' output(2,:)] )
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
%      From : 2001-DEC-14 20:10:14.196213 (TDB)
%      To   : 2001-DEC-14 21:35:50.318416 (TDB)
%
%
%   Example(2):
%
%      Find occultations of Titan by Saturn or of Saturn by
%      Titan as seen from the center of the Earth over the
%      last three months of 2008. Search for every type
%      of occultation.
%
%      Use light time corrections to model apparent positions of
%      Saturn and Titan. Stellar aberration corrections are not
%      specified because they don't affect occultation computations.
%
%      We select a step size of 15 minutes, which means we
%      ignore occultation events lasting less than 15 minutes,
%      if any exist.
%
%      MAXWIN  =  1000;
%      TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%      OCCTYP  = {'FULL', 'ANNULAR', 'PARTIAL', 'ANY' };
%
%
%      %
%      % Load kernels.
%      %
%      cspice_furnsh( 'standard.tm' );
%      cspice_furnsh( 'sat288.bsp' );
%
%
%      %
%      % Store the time bounds of our search interval in
%      % the cnfine confinement window.
%      %
%      et = cspice_str2et( { '2008 SEP 01 00:00:00 TDB', ...
%                            '2009 JAN 01 00:00:00 TDB'} );
%
%      cnfine = cspice_wninsd( et(1), et(2) );
%
%      %
%      % Select a 15-minute step. We'll ignore any occultations
%      % lasting less than 15 minutes.
%      %
%      step    = 900.;
%
%      %
%      % The observation location is the Earth.
%      %
%      obsrvr  = 'earth';
%      shape   = 'ellipsoid';
%      abcorr = 'lt';
%
%      for i=1:numel(OCCTYP)
%
%         %
%         % For each type, do a search for both transits of
%         % Titan across Saturn and occultations of Titan by
%         % Saturn.
%         %
%         for j=1:2
%
%            if isequal(j,1)
%
%               front  = 'TITAN';
%               fframe = 'IAU_TITAN';
%               back   = 'SATURN';
%               bframe = 'IAU_SATURN';
%
%            else
%
%               front  = 'SATURN';
%               fframe = 'IAU_SATURN';
%               back   = 'TITAN';
%               bframe = 'IAU_TITAN';
%
%            end
%
%            result = cspice_gfoclt( OCCTYP(i),                    ...
%                                    front, shape, fframe,         ...
%                                    back, shape, bframe,          ...
%                                    abcorr, obsrvr, step, cnfine, ...
%                                    MAXWIN);
%
%
%            fprintf( 'Condition      : %s\n',   char(OCCTYP(i)) )
%            fprintf( 'Occultation of : %s\n',   back  )
%            fprintf( 'by             : %s\n\n', front )
%
%
%            %
%            % List the beginning and ending times in each interval
%            % if result contains data.
%            %
%
%            count = numel(result)/2;
%
%            if isequal(count,0)
%
%               fprintf( 'Result window is empty.\n\n' )
%
%            else
%
%               for k=1:count
%
%                  [left, right] = cspice_wnfetd( result, k );
%
%                  output = cspice_timout( [left,right], TIMFMT );
%
%                  if( isequal( left, right) )
%
%                     disp( ['Event time: ' output(1,:)] )
%
%                  else
%
%                     disp( ['From : ' output(1,:)] )
%                     disp( ['To   : ' output(2,:)] )
%                     disp( ' ')
%
%                  end
%
%               end
%
%            end
%
%            %
%            % We've finished displaying the results of the
%            % current search.
%            %
%
%         end
%
%         %
%         % We've finished displaying the results of the
%         % searches using the current occultation type.
%         %
%
%      end
%
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Condition      : FULL
%      Occultation of : SATURN
%      by             : TITAN
%
%      Result window is empty.
%
%      Condition      : FULL
%      Occultation of : TITAN
%      by             : SATURN
%
%      From : 2008-OCT-27 22:08:01.627053 (TDB)
%      To   : 2008-OCT-28 01:05:03.375237 (TDB)
%
%      From : 2008-NOV-12 21:21:59.252263 (TDB)
%      To   : 2008-NOV-13 02:06:05.053051 (TDB)
%
%      From : 2008-NOV-28 20:49:02.402832 (TDB)
%      To   : 2008-NOV-29 02:13:58.986344 (TDB)
%
%      From : 2008-DEC-14 20:05:09.246177 (TDB)
%      To   : 2008-DEC-15 01:44:53.523002 (TDB)
%
%      From : 2008-DEC-30 19:00:56.577073 (TDB)
%      To   : 2008-DEC-31 00:42:43.222909 (TDB)
%
%      Condition      : ANNULAR
%      Occultation of : SATURN
%      by             : TITAN
%
%      From : 2008-OCT-19 21:29:20.599088 (TDB)
%      To   : 2008-OCT-19 22:53:34.518737 (TDB)
%
%      From : 2008-NOV-04 20:15:38.620369 (TDB)
%      To   : 2008-NOV-05 00:18:59.139979 (TDB)
%
%      From : 2008-NOV-20 19:38:59.647712 (TDB)
%      To   : 2008-NOV-21 00:35:26.725909 (TDB)
%
%      From : 2008-DEC-06 18:58:34.073269 (TDB)
%      To   : 2008-DEC-07 00:16:17.647040 (TDB)
%
%      From : 2008-DEC-22 18:02:46.288290 (TDB)
%      To   : 2008-DEC-22 23:26:52.712459 (TDB)
%
%      Condition      : ANNULAR
%      Occultation of : TITAN
%      by             : SATURN
%
%      Result window is empty.
%
%      Condition      : PARTIAL
%      Occultation of : SATURN
%      by             : TITAN
%
%      From : 2008-OCT-19 20:44:30.326772 (TDB)
%      To   : 2008-OCT-19 21:29:20.599088 (TDB)
%
%      From : 2008-OCT-19 22:53:34.518737 (TDB)
%      To   : 2008-OCT-19 23:38:26.250580 (TDB)
%
%      From : 2008-NOV-04 19:54:40.339331 (TDB)
%      To   : 2008-NOV-04 20:15:38.620369 (TDB)
%
%      From : 2008-NOV-05 00:18:59.139979 (TDB)
%      To   : 2008-NOV-05 00:39:58.612936 (TDB)
%
%      From : 2008-NOV-20 19:21:46.689523 (TDB)
%      To   : 2008-NOV-20 19:38:59.647712 (TDB)
%
%      From : 2008-NOV-21 00:35:26.725909 (TDB)
%      To   : 2008-NOV-21 00:52:40.604704 (TDB)
%
%      From : 2008-DEC-06 18:42:36.100544 (TDB)
%      To   : 2008-DEC-06 18:58:34.073269 (TDB)
%
%      From : 2008-DEC-07 00:16:17.647040 (TDB)
%      To   : 2008-DEC-07 00:32:16.324244 (TDB)
%
%      From : 2008-DEC-22 17:47:10.776723 (TDB)
%      To   : 2008-DEC-22 18:02:46.288290 (TDB)
%
%      From : 2008-DEC-22 23:26:52.712459 (TDB)
%      To   : 2008-DEC-22 23:42:28.850543 (TDB)
%
%      Condition      : PARTIAL
%      Occultation of : TITAN
%      by             : SATURN
%
%      From : 2008-OCT-27 21:37:16.970175 (TDB)
%      To   : 2008-OCT-27 22:08:01.627053 (TDB)
%
%      From : 2008-OCT-28 01:05:03.375237 (TDB)
%      To   : 2008-OCT-28 01:35:49.266507 (TDB)
%
%      From : 2008-NOV-12 21:01:47.105499 (TDB)
%      To   : 2008-NOV-12 21:21:59.252263 (TDB)
%
%      From : 2008-NOV-13 02:06:05.053051 (TDB)
%      To   : 2008-NOV-13 02:26:18.227358 (TDB)
%
%      From : 2008-NOV-28 20:31:28.522707 (TDB)
%      To   : 2008-NOV-28 20:49:02.402832 (TDB)
%
%      From : 2008-NOV-29 02:13:58.986344 (TDB)
%      To   : 2008-NOV-29 02:31:33.691598 (TDB)
%
%      From : 2008-DEC-14 19:48:27.094229 (TDB)
%      To   : 2008-DEC-14 20:05:09.246177 (TDB)
%
%      From : 2008-DEC-15 01:44:53.523002 (TDB)
%      To   : 2008-DEC-15 02:01:36.360243 (TDB)
%
%      From : 2008-DEC-30 18:44:23.485899 (TDB)
%      To   : 2008-DEC-30 19:00:56.577073 (TDB)
%
%      From : 2008-DEC-31 00:42:43.222909 (TDB)
%      To   : 2008-DEC-31 00:59:17.030569 (TDB)
%
%      Condition      : ANY
%      Occultation of : SATURN
%      by             : TITAN
%
%      From : 2008-OCT-19 20:44:30.326772 (TDB)
%      To   : 2008-OCT-19 23:38:26.250580 (TDB)
%
%      From : 2008-NOV-04 19:54:40.339331 (TDB)
%      To   : 2008-NOV-05 00:39:58.612936 (TDB)
%
%      From : 2008-NOV-20 19:21:46.689523 (TDB)
%      To   : 2008-NOV-21 00:52:40.604704 (TDB)
%
%      From : 2008-DEC-06 18:42:36.100544 (TDB)
%      To   : 2008-DEC-07 00:32:16.324244 (TDB)
%
%      From : 2008-DEC-22 17:47:10.776723 (TDB)
%      To   : 2008-DEC-22 23:42:28.850543 (TDB)
%
%      Condition      : ANY
%      Occultation of : TITAN
%      by             : SATURN
%
%      From : 2008-OCT-27 21:37:16.970175 (TDB)
%      To   : 2008-OCT-28 01:35:49.266507 (TDB)
%
%      From : 2008-NOV-12 21:01:47.105499 (TDB)
%      To   : 2008-NOV-13 02:26:18.227358 (TDB)
%
%      From : 2008-NOV-28 20:31:28.522707 (TDB)
%      To   : 2008-NOV-29 02:31:33.691598 (TDB)
%
%      From : 2008-DEC-14 19:48:27.094229 (TDB)
%      To   : 2008-DEC-15 02:01:36.360243 (TDB)
%
%      From : 2008-DEC-30 18:44:23.485899 (TDB)
%      To   : 2008-DEC-31 00:59:17.030569 (TDB)
%
%-Particulars
%
%   This routine provides a simple interface for conducting searches for
%   occultation events.
%
%   This routine determines a set of one or more time intervals
%   within the confinement window when a specified type of
%   occultation occurs. The resulting set of intervals is returned as
%   a SPICE window.
%
%   Below we discuss in greater detail aspects of this routine's
%   solution process that are relevant to correct and efficient
%   use of this routine in user applications.
%
%   The Search Process
%   ==================
%
%   The search for occultations is treated as a search for state
%   transitions: times are sought when the state of the `back' body
%   changes from "not occulted" to "occulted" or vice versa.
%
%   Step Size
%   =========
%
%   Each interval of the confinement window is searched as follows:
%   first, the input step size is used to determine the time
%   separation at which the occultation state will be sampled.
%   Starting at the left endpoint of an interval, samples will be
%   taken at each step. If a state change is detected, a root has
%   been bracketed; at that point, the "root"--the time at which the
%   state change occurs---is found by a refinement process, for
%   example, by a binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the occultation state is constant:
%   the step size should be shorter than the shortest occultation
%   duration and the shortest period between occultations, within
%   the confinement window.
%
%   Having some knowledge of the relative geometry of the targets and
%   observer can be a valuable aid in picking a reasonable step size.
%   In general, the user can compensate for lack of such knowledge by
%   picking a very short step size; the cost is increased computation
%   time.
%
%   Note that the step size is not related to the precision with which
%   the endpoints of the intervals of the result window are computed.
%   That precision level is controlled by the convergence tolerance.
%
%   Convergence Tolerance
%   =====================
%
%   As described above, the root-finding process used by this routine
%   involves first bracketing roots and then using a search process
%   to locate them. "Roots" are both times when local extrema are
%   attained and times when the distance function is equal to a
%   reference value. All endpoints of the intervals comprising the
%   result window are either endpoints of intervals of the
%   confinement window or roots.
%
%   Once a root has been bracketed, a refinement process is used to
%   narrow down the time interval within which the root must lie.
%   This refinement process terminates when the location of the root
%   has been determined to within an error margin called the
%   "convergence tolerance." The convergence tolerance used by this
%   routine is set by the parameter SPICE_GF_CNVTOL.
%
%   The value of SPICE_GF_CNVTOL is set to a "tight" value so that the
%   tolerance doesn't become the limiting factor in the accuracy of
%   solutions found by this routine. In general the accuracy of input
%   data will be the limiting factor.
%
%   The user may change the convergence tolerance from the default
%   SPICE_GF_CNVTOL value by calling the routine cspice_gfstol, e.g.
%
%      cspice_gfstol( tolerance value in seconds )
%
%   Call cspice_gfstol prior to calling this routine. All subsequent
%   searches will use the updated tolerance value.
%
%   Setting the tolerance tighter than SPICE_GF_CNVTOL is unlikely to be
%   useful, since the results are unlikely to be more accurate.
%   Making the tolerance looser will speed up searches somewhat,
%   since a few convergence steps will be omitted. However, in most
%   cases, the step size is likely to have a much greater affect on
%   processing time than would the convergence tolerance.
%
%   The Confinement Window
%   ======================
%
%   The simplest use of the confinement window is to specify a time
%   interval within which a solution is sought. However, the
%   confinement window can, in some cases, be used to make searches
%   more efficient. Sometimes it's possible to do an efficient search
%   to reduce the size of the time period over which a relatively
%   slow search of interest must be performed.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine gfoclt_c.
%
%   MICE.REQ
%   GF.REQ
%   SPK.REQ
%   CK.REQ
%   TIME.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.1.0, 12-MAY-2012, EDW (JPL)
%
%      Renamed the argument 'size' to 'room'. "size" is a Matlab function
%      name and it's seriously dumb to use a function name word as an argument
%      name.
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%      Header updated to describe use of cspice_gfstol.
%
%   -Mice Version 1.0.0, 15-APR-2009, EDW (JPL)
%
%-Index_Entries
%
%   GF occultation search
%
%-&

function [result] = cspice_gfoclt( occtyp, front,  fshape, fframe, ...
                                   back,   bshape, bframe, abcorr, ...
                                   obsrvr, step,   cnfine, room)

   switch nargin

      case 12

         occtyp  = zzmice_str(occtyp);
         front   = zzmice_str(front);
         fshape  = zzmice_str(fshape);
         fframe  = zzmice_str(fframe);
         back    = zzmice_str(back);
         bshape  = zzmice_str(bshape);
         bframe  = zzmice_str(bframe);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         step    = zzmice_dp(step);
         cnfine  = zzmice_win(cnfine);
         room    = zzmice_int(room,    [1, int32(inf)/2] );

      otherwise

         error ( [ 'Usage: [result] = cspice_gfoclt( `occtyp`, '          ...
                                           '`front`, `fshape`, `fframe`, '...
                                           '`back`, `bshape`, `bframe`, ' ...
                                           '`abcorr`, `obsrvr`, step, '   ...
                                           'cnfine, room)' ] )

   end

   %
   % Call the GF routine.  Add to 'cnfine' the 6x1 space needed for
   % the control segment.
   %
   try

      [result] = mice('gfoclt_c',  occtyp, front,  fshape, fframe, ...
                                   back,   bshape, bframe,         ...
                                   abcorr, obsrvr,  step,          ...
                                   [zeros(6,1); cnfine], room);
   catch
      rethrow(lasterror)
   end




