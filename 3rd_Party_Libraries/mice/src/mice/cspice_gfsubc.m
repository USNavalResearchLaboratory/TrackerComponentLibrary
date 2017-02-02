%-Abstract
%
%   CSPICE_GFSUBC determines the time intervals for which a coordinate
%   of an subpoint position vector satisfies a numerical constraint.
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
%     target   the string scalar naming the target body.  Optionally,
%               you may supply the integer ID code for the object as an
%               integer string.  For example both 'MOON' and '301'
%               are legitimate strings that indicate the moon is the
%               target body.
%
%               On calling cspice_gfsubc, the kernel pool must contain the
%               radii data corresponding to 'target'.
%
%      fixref   the string scalar naming the body-fixed, body-centered
%               reference frame associated with the target body 'target'.
%
%               The SPICE frame subsystem must recognize the 'fixref' name.
%
%      method   the string scalar naming the method to use for the subpoint
%               calculation. The accepted values for method:
%
%                 'Near point: ellipsoid'   The sub-observer point
%                                           computation uses a
%                                           triaxial ellipsoid to
%                                           model the surface of the
%                                           target body. The
%                                           sub-observer point is
%                                           defined as the nearest
%                                           point on the target
%                                           relative to the
%                                           observer.
%
%                 'Intercept: ellipsoid'    The sub-observer point
%                                           computation uses a
%                                           triaxial ellipsoid to
%                                           model the surface of the
%                                           target body. The
%                                           sub-observer point is
%                                           defined as the target
%                                           surface intercept of the
%                                           line containing the
%                                           observer and the
%                                           target's center.
%

%               The 'method' string lacks sensitivity to case, and to embedded,
%               leading and trailing blanks.
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
%      crdsys   the string scalar naming the coordinate system for which the
%               coordinate of interest is a member.
%
%      coord    the string scalar naming the coordinate of interest in 'crdsys'.
%
%               The supported coordinate systems and coordinate names are:
%
%               Coordinate System (crdsys)    Coordinates (coord)      Range
%
%                 'RECTANGULAR'                  'X'
%                                                'Y'
%                                                'Z'
%
%                 'LATITUDINAL'                  'RADIUS'
%                                                'LONGITUDE'        (-Pi,Pi]
%                                                'LATITUDE'         [-Pi/2,Pi/2]
%
%                 'RA/DEC'                       'RANGE'
%                                                'RIGHT ASCENSION'  [0,2Pi)
%                                                'DECLINATION'      [-Pi/2,Pi/2]
%
%                 'SPHERICAL'                    'RADIUS'
%                                                'COLATITUDE'       [0,Pi]
%                                                'LONGITUDE'        (-Pi,Pi]
%
%                 'CYLINDRICAL'                  'RADIUS'
%                                                'LONGITUDE'        [0,2Pi)
%                                                'Z'
%
%                 'GEODETIC'                     'LONGITUDE'        (-Pi,Pi]
%                                                'LATITUDE'         [-Pi/2,Pi/2]
%                                                'ALTITUDE'
%
%                 'PLANETOGRAPHIC'               'LONGITUDE'        [0,2Pi)
%                                                'LATITUDE'         [-Pi/2,Pi/2]
%                                                'ALTITUDE'
%
%                 The ALTITUDE coordinates have a constant value
%                 of zero +/- roundoff for ellipsoid targets.
%
%                 Limit searches for coordinate events in the GEODETIC and
%                 PLANETOGRAPHIC coordinate systems to 'target' bodies with
%                 axial symmetry in the equatorial plane, i.e. equality
%                 of the body X and Y radii (oblate or prolate spheroids).
%
%      relate   the string scalar or character describing the relational
%               operator used to define a constraint on the selected coordinate
%               of the subpoint vector. The result window found by this routine
%               indicates the time intervals where the constraint is satisfied.
%               Supported values of relate and corresponding meanings are
%               shown below:
%
%                  '>'       Separation is greater than the reference
%                            value refval.
%
%                  '='       Separation is equal to the reference
%                            value refval.
%
%                  '<'      Separation is less than the reference
%                            value refval.
%
%                  'ABSMAX'  Separation is at an absolute maximum.
%
%                  'ABSMIN'  Separation is at an absolute  minimum.
%
%                  'LOCMAX'  Separation is at a local maximum.
%
%                  'LOCMIN'  Separation is at a local minimum.
%
%               The caller may indicate that the region of interest
%               is the set of time intervals where the quantity is
%               within a specified measure of an absolute extremum.
%               The argument 'adjust' (described below) is used to
%               specify this measure.
%
%               Local extrema are considered to exist only in the
%               interiors of the intervals comprising the confinement
%               window:  a local extremum cannot exist at a boundary
%               point of the confinement window.
%
%               The 'relate' string lacks sensitivity to case, and to leading
%               and trailing blanks.
%
%      refval   the double precision scalar reference value used together
%               with relate argument to define an equality or inequality to
%               satisfy by the selected coordinate of the subpoint vector.
%               See the discussion of relate above for further
%               information.
%
%               The units of 'refval' correspond to the type as defined
%               by 'coord', radians for angular measures, kilometers for
%               distance measures.
%
%      adjust   a double precision scalar value used to modify searches for
%               absolute extrema: when relate is set to ABSMAX or ABSMIN and
%               adjust is set to a positive value, cspice_gfsubc finds times
%               when the subpoint vector coordinate is within 'adjust'
%               radians/kilometers of the specified extreme value.
%
%               For relate set to ABSMAX, the result window contains
%               time intervals when the subpoint vector coordinate has
%               values between ABSMAX - 'adjust' and ABSMAX.
%
%               For relate set to ABSMIN, the result window contains
%               time intervals when the position vector coordinate has
%               values between ABSMIN and ABSMIN + 'adjust'.
%
%               'adjust' is not used for searches for local extrema,
%               equality or inequality conditions.
%
%      step     the double precision time step size to use in the search.
%
%               'step' must be short enough to for a search using this step
%               size to locate the time intervals where coordinate
%               function of the subpoint vector is monotone increasing or
%               decreasing. However, 'step must not be *too* short, or
%               the search will take an unreasonable amount of time.
%
%               For coordinates other than LONGITUDE and RIGHT ASCENSION,
%               the step size must be shorter than the shortest interval,
%               within the confinement window, over which the coordinate
%               is monotone increasing or decreasing.
%
%               For LONGITUDE and RIGHT ASCENSION, the step size must
%               be shorter than the shortest interval, within the
%               confinement window, over which either the sin or cos
%               of the coordinate is monotone increasing or decreasing.
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
%      result = cspice_gfsubc( target, fixref, method, abcorr, ...
%                              obsrvr, crdsys, coord, relate,  ...
%                              refval, adjust, step,  nintvls, ...
%                              cnfine )
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
%      Find the time during 2007 for which the subpoint position vector
%      of the sun on earth in the IAU_EARTH frame lies within a geodetic
%      latitude-longitude "box" defined as
%
%         16 degrees <= latitude  <= 17 degrees
%         85 degrees <= longitude <= 86 degrees
%
%      This problem requires four searches, each search on one of the
%      box restrictions. The user needs also realize the temporal
%      behavior of latitude greatly differs from that of the longitude. The
%      sub-observer point latitude varies between approximately 23.44 degrees
%      and -23.44 degrees during the year. The sub-observer point longitude
%      varies between -180 degrees and 180 degrees in one day.
%
%      MAXWIN  =  1000;
%      TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%      R2D     = cspice_dpr;
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
%      % The latitude varies relatively slowly (46 degrees) during the
%      % year. The extrema occur approximately every six months.
%      % Search using a step size less than half that value (180 days).
%      % For this example use ninety days (in units of seconds).
%      %
%      step   = cspice_spd*90.;
%
%      %
%      % Perform four searches to determine the times when the latitude-
%      % longitude box restriction conditions apply to the subpoint vector.
%      %
%      % Use geodetic coordinates.
%      %
%      adjust = 0.;
%      target = 'EARTH';
%      obsrvr = 'SUN';
%      method = 'Near point: ellipsoid';
%      fixref = 'IAU_EARTH';
%      crdsys = 'GEODETIC';
%      nintvls= MAXWIN;
%      abcorr = 'NONE';
%
%      %
%      %  Perform the searches such that the result window of a search
%      % serves as the confinement window of the subsequent search.
%      %
%
%      %
%      % Since the latitude coordinate varies slowly and is well behaved
%      % over the time of the confinement window, search first for the
%      % windows satisfying the latitude requirements, then use that result
%      % as confinement for the longitude search.
%      %
%      coord  = 'LATITUDE';
%      refval = 16. * cspice_rpd;
%      relate = '>';
%
%      %
%      % Perform this search using the geometric position
%      % of the bodies; set the aberration correction to 'NONE'.
%      %
%
%      result1 = cspice_gfsubc( target, fixref, method, abcorr, obsrvr, ...
%                              crdsys, coord, relate, refval, ...
%                              adjust, step, nintvls, cnfine );
%
%      refval = 17. * cspice_rpd;
%      relate = '<';
%
%      result2 = cspice_gfsubc( target, fixref, method, abcorr, obsrvr, ...
%                               crdsys, coord, relate, refval,  ...
%                               adjust, step, nintvls, result1 );
%
%      %
%      % Now the longitude search.
%      %
%      coord  = 'LONGITUDE';
%
%      %
%      % Reset the step size to something appropriate for the 360
%      % degrees in 24 hours domain. The longitude shows near
%      % linear behavior so use a step size less than half the period
%      % of twelve hours. Ten hours will suffice in this case.
%      %
%      step   = cspice_spd * (10./24.);
%
%      refval = 85. * cspice_rpd;
%      relate = '>';
%
%      result3 = cspice_gfsubc( target, fixref, method, abcorr, obsrvr, ...
%                               crdsys, coord, relate, refval,  ...
%                               adjust, step, nintvls, result2 );
%
%      %
%      % Contract the endpoints of each window to account
%      % for possible round-off error at the -180/180 degree branch.
%      %
%      % A contraction value of a millisecond should eliminate
%      % any round-off caused branch crossing.
%      %
%      result3 = cspice_wncond( 1e-3, 1e-3, result3 );
%
%      refval = 86. * cspice_rpd;
%      relate = '<';
%
%      result4 = cspice_gfsubc( target, fixref, method, abcorr, obsrvr, ...
%                               crdsys, coord, relate, refval,  ...
%                               adjust, step, nintvls, result3 );
%
%      result = result4;
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
%            [lspoint, trgepc, srfvec] = cspice_subpnt(method, target, left, ...
%                                                      fixref, abcorr, obsrvr);
%            [lrad, llong, llat ] = cspice_reclat( lspoint );
%
%            [rspoint, trgepc, srfvec] = cspice_subpnt(method, target,  ...
%                                                right, fixref, abcorr, obsrvr);
%            [rrad, rlong, rlat ] = cspice_reclat( rspoint );
%
%            fprintf( 'From : %s   %f   %f\n', output(1,:), ...
%                                              llat*R2D,    ...
%                                              llong*R2D )
%            fprintf( 'To   : %s   %f   %f\n', output(2,:), ...
%                                              rlat*R2D,    ...
%                                              rlong*R2D )
%            fprintf( '\n' );
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
%      From : 2007-MAY-05 06:14:04.637735 (TDB)   16.054356   86.000000
%      To   : 2007-MAY-05 06:18:04.621907 (TDB)   16.055148   85.000004
%
%      From : 2007-MAY-06 06:13:59.583483 (TDB)   16.337147   86.000000
%      To   : 2007-MAY-06 06:17:59.569240 (TDB)   16.337927   85.000004
%
%      From : 2007-MAY-07 06:13:55.102940 (TDB)   16.615444   86.000000
%      To   : 2007-MAY-07 06:17:55.090299 (TDB)   16.616210   85.000004
%
%      From : 2007-MAY-08 06:13:51.202604 (TDB)   16.889163   86.000000
%      To   : 2007-MAY-08 06:17:51.191583 (TDB)   16.889916   85.000004
%
%      From : 2007-AUG-06 06:23:17.282927 (TDB)   16.680717   86.000000
%      To   : 2007-AUG-06 06:27:17.264010 (TDB)   16.679962   85.000004
%
%      From : 2007-AUG-07 06:23:10.545441 (TDB)   16.406411   86.000000
%      To   : 2007-AUG-07 06:27:10.524925 (TDB)   16.405643   85.000004
%
%      From : 2007-AUG-08 06:23:03.233996 (TDB)   16.127678   86.000000
%      To   : 2007-AUG-08 06:27:03.211889 (TDB)   16.126897   85.000004
%
%-Particulars
%
%   This routine provides a simple interface for conducting searches
%   for subpoint surface vector coordinate value events.
%
%   This routine determines a set of one or more time intervals
%   within the confinement window when the selected coordinate of
%   the subpoint surface vector satisfies a caller-specified
%   constraint. The resulting set of intervals is returned as a SPICE
%   window.
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
%   coordinate function is monotone increasing and monotone
%   decreasing. Each of these time periods is represented by a SPICE
%   window. Having found these windows, all of the coordinate
%   function's local extrema within the confinement window are known.
%   Absolute extrema then can be found very easily.
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
%   change of coordinate will be sampled. Starting at
%   the left endpoint of an interval, samples will be taken at each
%   step. If a change of sign is found, a root has been bracketed; at
%   that point, the time at which the time derivative of the coordinate
%   is zero can be found by a refinement process, for example,
%   using a binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the coordinate function is monotone:
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
%   Practical use of the coordinate search capability would likely
%   consist of searches over multiple coordinate constraints to find
%   time intervals that satisfies the constraints. An effective
%   technique to accomplish such a search is to use the result
%   window from one search as the confinement window of the next.
%
%   Longitude and Right Ascension
%   =============================
%
%   The cyclic nature of the longitude and right ascension coordinates
%   produces branch cuts at +/- 180 degrees longitude and 0-360
%   longitude. Round-off error may cause solutions near these branches
%   to cross the branch. Use of the SPICE routine wncond_c will contract
%   solution windows by some epsilon, reducing the measure of the
%   windows and eliminating the branch crossing. A one millisecond
%   contraction will in most cases eliminate numerical round-off caused
%   branch crossings.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine gfsubc_c.
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
%      Minor header edit "cosin" -> "cos".
%
%      Edit to Example description, replaced "intercept" with
%      "sub-observer point."
%
%      Header updated to describe use of cspice_gfstol.
%
%   -Mice Version 1.0.0, 15-APR-2009, EDW (JPL)
%
%-Index_Entries
%
%   GF subpoint coordinate search
%
%-&

function [result] = cspice_gfsubc( target, fixref, method, abcorr, ...
                                   obsrvr, crdsys, coord, relate,  ...
                                   refval, adjust, step,  nintvls, ...
                                   cnfine )

   switch nargin

      case 13

         target  = zzmice_str(target);
         fixref  = zzmice_str(fixref);
         method  = zzmice_str(method);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         crdsys  = zzmice_str(crdsys);
         coord   = zzmice_str(coord);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);

      otherwise

         error ( [ 'Usage: [result] = cspice_gfsubc( `target`, `fixref`, '  ...
                                 '`method`, `abcorr`, `obsrvr`, `crdsys`, ' ...
                                 '`coord`, `relate`, refval, adjust, '      ...
                                 'step, nintvls, cnfine )' ] )

   end

   %
   % Call the GF routine, add to 'cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gfsubc_c', target, fixref, method, abcorr, ...
                                  obsrvr, crdsys, coord, relate,  ...
                                  refval, adjust, step,  nintvls, ...
                                  [zeros(6,1); cnfine] );
   catch
      rethrow(lasterror)
   end




