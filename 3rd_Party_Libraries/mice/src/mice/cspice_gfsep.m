%-Abstract
%
%   CSPICE_GFSEP determines the time intervals when the angular separation
%   between the position vectors of two target bodies relative to an
%   observer satisfies a numerical relationship.
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
%      targ1    the scalar string naming the first body of interest. You can
%               also supply the integer ID code for the object as an
%               integer string.  For example both 'MOON' and '301'
%               are legitimate strings that indicate the moon is the
%               target body.
%
%      shape1   the scalar string naming the geometric model used to represent
%               the shape of the 'targ1' body. Models supported by this routine:
%
%                  'SPHERE'        Treat the body as a sphere with radius
%                                  equal to the maximum value of
%                                  BODYnnn_RADII
%
%                  'POINT'         Treat the body as a point;
%                                  radius has value zero.
%
%                  The 'shape1' string lacks sensitivity to case, leading
%                  and trailing blanks.
%
%      frame1   the scalar string naming the body-fixed reference frame
%               corresponding to 'targ1'. cspice_gfsep does not currently use
%               this argument's value, its use is reserved for future
%               shape models. The value 'NULL' will suffice for
%               "POINT" and "SPHERE" shaped bodies.
%
%      targ2    the scalar string naming the second body of interest. You can
%               also supply the integer ID code for the object as an
%               integer string.  For example both 'MOON' and '301'
%               are legitimate strings that indicate the moon is the
%               target body.
%
%      shape2   the scalar string naming the geometric model used to represent
%               the shape of the 'targ2.' Models supported by this routine:
%
%                 'SPHERE'        Treat the body as a sphere with radius
%                                 equal to the maximum value of
%                                 BODYnnn_RADII
%
%                 'POINT'         Treat the body as a single point;
%                                 radius has value zero.
%
%               The 'shape2' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      frame2   the scalar string naming the body-fixed reference frame
%               corresponding to 'targ2'. cspice_gfsep does not currently use
%               this argument's value, its use is reserved for future
%               shape models. The value 'NULL' will suffice for
%               "POINT" and "SPHERE" shaped bodies.
%
%      abcorr   the string scalar indicating the aberration corrections to apply
%               to the state evaluations to account for one-way light time and
%               stellar aberration.
%
%               This routine accepts the same aberration corrections as does
%               the routine spkezr_c. See the header of spkezr_c for a
%               detailed description of the aberration correction options.
%               For convenience, the options are listed below:
%
%               This routine accepts the same aberration corrections as does
%               the SPICE routine SPKEZR. See the header of SPKEZR for a
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
%      obsrvr   the string scalar naming the observing body. Optionally, you
%               may supply the ID code of the object as an integer
%               string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the
%               observer is Earth.
%
%      relate   the string identifying the relational operator used to
%               define a constraint on the angular separation. The result
%               window found by this routine indicates the time intervals
%               where the constraint is satisfied. Supported values of
%               relate and corresponding meanings are shown below:
%
%                  '>'      Separation is greater than the reference
%                           value 'refval'.
%
%                  '='      Separation is equal to the reference
%                           value 'refval'.
%
%                  '<'      Separation is less than the reference
%                           value 'refval'.
%
%                 'ABSMAX'  Separation is at an absolute maximum.
%
%                 'ABSMIN'  Separation is at an absolute  minimum.
%
%                 'LOCMAX'  Separation is at a local maximum.
%
%                 'LOCMIN'  Separation is at a local minimum.
%
%               The caller may indicate that the region of interest
%               is the set of time intervals where the quantity is
%               within a specified angular separation of an absolute extremum.
%               The argument adjust (described below) is used to
%               specify this angular separation.
%
%               Local extrema are considered to exist only in the
%               interiors of the intervals comprising the confinement
%               window:  a local extremum cannot exist at a boundary
%               point of the confinement window.
%
%               Negative Angular Separation
%
%                  For those searches using a SPHERE shape identifier for
%                  either target body, the angular separation function
%                  returns a negative value when the bodies overlap (occult).
%
%               The 'relate' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      refval   the double precision reference value used together with
%               'relate' argument to define an equality or inequality to be
%               satisfied by the angular separation between the specified target
%               and observer. See the discussion of relate above for
%               further information.
%
%               The units of 'refval' are radians.
%
%      adjust   a double precision value used to modify searches for
%               absolute extrema: when 'relate' is set to ABSMAX or ABSMIN and
%               adjust is set to a positive value, cspice_gfsep finds times when
%               the angular separation between the bodies is within adjust
%               radians of the specified extreme value.
%
%               For relate set to ABSMAX, the result window contains
%               time intervals when the angular separation has
%               values between ABSMAX - 'adjust' and ABSMAX.
%
%               For relate set to ABSMIN, the result window contains
%               time intervals when the angular separation has
%               values between ABSMIN and ABSMIN + 'adjust'.
%
%               'adjust' is not used for searches for local extrema,
%               equality or inequality conditions.
%
%      step     the double precision time step size to use in the search.
%
%               'step' must be short enough to for a search using this step
%               size to locate the time intervals where the
%               specified angular separation function is monotone
%               increasing or decreasing. However, 'step' must not be
%               *too* short, or the search will take an unreasonable
%               amount of time.
%
%               The choice of step affects the completeness but not
%               the precision of solutions found by this routine; the
%               precision is controlled by the convergence tolerance.
%               See the discussion of the parameter SPICE_GF_CNVTOL for
%               details.
%
%               'step' has units of TDB seconds.
%
%      nintvls  an integer scalar value specifying the number of intervals in
%               the internal workspace array used by this routine. 'nintvls'
%               should be at least as large as the number of intervals
%               within the search region on which the specified observer-target
%               vector coordinate function is monotone increasing or decreasing.
%               It does no harm to pick a value of 'nintvls' larger than the
%               minimum required to execute the specified search, but if chosen
%               too small, the search will fail.
%
%      cnfine   a double precision SPICE window that confines the time
%               period over which the specified search is conducted.
%               cnfine may consist of a single interval or a collection
%               of intervals.
%
%               In some cases the confinement window can be used to
%               greatly reduce the time period that must be searched
%               for the desired solution. See the Particulars section
%               below for further discussion.
%
%               See the Examples section below for a code example
%               that shows how to create a confinement window.
%
%   the call:
%
%      result = cspice_gfsep( targ1,  shape1, frame1,          ...
%                             targ2,  shape2, frame2,          ...
%                             abcorr, obsrvr, relate,  refval, ...
%                             adjust, step,   nintvls, cnfine )
%
%   returns:
%
%      result   the SPICE window of intervals, contained within the
%               confinement window 'cnfine', on which the specified
%               constraint is satisfied.
%
%               If the search is for local extrema, or for absolute
%               extrema with adjust set to zero, then normally each
%               interval of result will be a singleton: the left and
%               right endpoints of each interval will be identical.
%
%               If no times within the confinement window satisfy the
%               constraint, 'result' will return with cardinality zero.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Determine the times of local maxima of the angular separation
%      between the moon and sun as observed from earth from
%      Jan 1, 2007 to Jan 1 2008.
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
%      et = cspice_str2et( { '2007 JAN 01', '2008 JAN 01'} );
%
%      cnfine = cspice_wninsd( et(1), et(2) );
%
%      %
%      % Search using a step size of 6 days (in units of seconds).
%      %
%      step   = 6.*cspice_spd;
%      adjust = 0.;
%      refval = 0;
%
%      targ1  = 'MOON';
%      shape1 = 'SPHERE';
%      frame1 = 'NULL';
%      targ2  = 'SUN';
%      shape2 = 'SPHERE';
%      frame2 = 'NULL';
%      abcorr = 'NONE';
%      relate = 'LOCMAX';
%      obsrvr = 'EARTH';
%      nintvls = MAXWIN;
%
%      result = cspice_gfsep( targ1,  shape1, frame1, ...
%                             targ2,  shape2, frame2, ...
%                             abcorr, obsrvr, relate, ...
%                             refval, adjust, step,   ...
%                             nintvls, cnfine );
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
%      Event time: 2007-JAN-03 14:20:24.618884 (TDB)
%      Event time: 2007-FEB-02 06:16:24.101655 (TDB)
%      Event time: 2007-MAR-03 23:22:41.994289 (TDB)
%      Event time: 2007-APR-02 16:49:16.134481 (TDB)
%      Event time: 2007-MAY-02 09:41:43.829169 (TDB)
%      Event time: 2007-JUN-01 01:03:44.527040 (TDB)
%      Event time: 2007-JUN-30 14:15:26.576639 (TDB)
%      Event time: 2007-JUL-30 01:14:49.002265 (TDB)
%      Event time: 2007-AUG-28 10:39:01.390508 (TDB)
%      Event time: 2007-SEP-26 19:25:51.512445 (TDB)
%      Event time: 2007-OCT-26 04:30:56.628530 (TDB)
%      Event time: 2007-NOV-24 14:31:04.334590 (TDB)
%      Event time: 2007-DEC-24 01:40:12.238389 (TDB)
%
%-Particulars
%
%   This routine provides a simple interface for conducting searches
%   for angular separation events.
%
%   This routine determines a set of one or more time intervals
%   within the confinement window for which the angular separation
%   between the two bodies satisfies some defined relationship.
%   The resulting set of intervals is returned as a SPICE window.
%
%   Below we discuss in greater detail aspects of this routine's
%   solution process that are relevant to correct and efficient
%   use of this routine in user applications.
%
%   The Search Process
%   ==================
%
%   Regardless of the type of constraint selected by the caller, this
%   routine starts the search for solutions by determining the time
%   periods, within the confinement window, over which the specified
%   angular separation function is monotone increasing and monotone
%   decreasing. Each of these time periods is represented by a SPICE window.
%   Having found these windows, all of the angular separation function's
%   local extrema within the confinement window are known. Absolute extrema
%   then can be found very easily.
%
%   Within any interval of these "monotone" windows, there will be at
%   most one solution of any equality constraint. Since the boundary
%   of the solution set for any inequality constraint is contained in
%   the union of
%
%      - the set of points where an equality constraint is met
%      - the boundary points of the confinement window
%
%   the solutions of both equality and inequality constraints can be
%   found easily once the monotone windows have been found.
%
%   Step Size
%   =========
%
%   The monotone windows (described above) are found using a two-step
%   search process. Each interval of the confinement window is
%   searched as follows: first, the input step size is used to
%   determine the time separation at which the sign of the rate of
%   change of angular separation (angular separation rate) will be
%   sampled. Starting at the left endpoint of an interval, samples
%   will be taken at each step. If a change of sign is found, a
%   root has been bracketed; at that point, the time at which the
%   angular separation rate is zero can be found by a refinement
%   process, for example, using a binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the distance function is monotone:
%   the step size should be shorter than the shortest of these
%   intervals (within the confinement window).
%
%   The optimal step size is *not* necessarily related to the lengths
%   of the intervals comprising the result window. For example, if
%   the shortest monotone interval has length 10 days, and if the
%   shortest result window interval has length 5 minutes, a step size
%   of 9.9 days is still adequate to find all of the intervals in the
%   result window. In situations like this, the technique of using
%   monotone windows yields a dramatic efficiency improvement over a
%   state-based search that simply tests at each step whether the
%   specified constraint is satisfied. The latter type of search can
%   miss solution intervals if the step size is longer than the
%   shortest solution interval.
%
%   Having some knowledge of the relative geometry of the target and
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
%   Elongation
%   ===========================
%
%   The angular separation of two targets as seen from an observer
%   where one of those targets is the sun is known as elongation.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine gfsep_c.
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
%   -Mice Version 1.0.1, 05-SEP-2012, EDW (JPL)
%
%      Edit to comments to correct search description.
%
%      Header updated to describe use of cspice_gfstol.
%
%    -Mice Version 1.0.1, 29-DEC-2009, EDW (JPL)
%
%      Edited argument descriptions. Removed mention of "ELLIPSOID"
%      shape from 'shape1' and 'shape2' as that option is not yet
%      implemented.
%
%   -Mice Version 1.0.0, 15-APR-2009, NJB EDW (JPL)
%
%-Index_Entries
%
%   GF angular separation search
%
%-&

function [result] = cspice_gfsep( targ1, shape1, frame1,          ...
                                  targ2, shape2, frame2,          ...
                                  abcorr, obsrvr, relate, refval, ...
                                  adjust, step, nintvls, cnfine )

   switch nargin

      case 14

         targ1   = zzmice_str(targ1);
         shape1  = zzmice_str(shape1);
         frame1  = zzmice_str(frame1);
         targ2   = zzmice_str(targ2);
         shape2  = zzmice_str(shape2);
         frame2  = zzmice_str(frame2);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);

      otherwise

         error ( [ 'Usage: [result] = cspice_gfsep( `targ1`, `shape1`, ' ...
                               '`frame1`, `targ2`, `shape2`, `frame2`, ' ...
                               '`abcorr`, `obsrvr`, `relate`, refval, '  ...
                               'adjust, step, nintvls, cnfine )' ] )

   end

   %
   % Call the GF routine, add to 'cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gfsep_c',  targ1, shape1, frame1,          ...
                                  targ2, shape2, frame2,          ...
                                  abcorr, obsrvr, relate, refval, ...
                                  adjust, step, nintvls,          ...
                                  [zeros(6,1); cnfine] );
   catch
      rethrow(lasterror)
   end




