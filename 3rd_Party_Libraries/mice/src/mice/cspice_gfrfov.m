%-Abstract
%
%   CSPICE_GFRFOV determine time intervals when a specified ray intersects the
%   space bounded by the field-of-view (FOV) of a specified instrument.
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
%      SPICE_GF_MAXVRT
%
%               is the maximum number of vertices that may be used
%               to define the boundary of the specified instrument's
%               field of view.
%
%      SPICE_GF_MARGIN
%
%               is a small positive number used to constrain the
%               orientation of the boundary vectors of polygonal
%               FOVs. Such FOVs must satisfy the following constraints:
%
%                  1) The boundary vectors must be contained within
%                  a right circular cone of angular radius less
%                  than (pi/2) - SPICE_GF_MARGIN radians; in other
%                  words, there must be a vector A such that all
%                  boundary vectors have angular separation from
%                  A of less than (pi/2)-SPICE_GF_MARGIN radians.
%
%                  2) There must be a pair of boundary vectors U, V
%                  such that all other boundary vectors lie in
%                  the same half space bounded by the plane
%                  containing U and V. Furthermore, all other
%                  boundary vectors must have orthogonal
%                  projections onto a plane normal to this plane
%                  such that the projections have angular
%                  separation of at least 2*SPICE_GF_MARGIN radians
%                  from the plane spanned by U and V.
%
%      Arguments-
%
%      inst     the string naming the instrument, such as a
%               spacecraft-mounted framing camera, the field of view
%               (FOV) of which is to be used for an target intersection
%               search: the direction from the observer to a target
%               is represented by a ray, and times when the specified
%               ray intersects the region of space bounded by the FOV
%               are sought.
%
%               The position of the instrument designated by 'inst' is
%               considered to coincide with that of the ephemeris
%               object designated by the input argument 'obsrvr' (see
%               description below).
%
%               'inst' must have a corresponding NAIF ID and a frame
%               defined, as is normally done in a frame kernel. It
%               must also have an associated reference frame and a FOV
%               shape, boresight and boundary vertices (or reference
%               vector and reference angles) defined, as is usually
%               done in an instrument kernel.
%
%               See the header of the Mice routine cspice_getfov for a
%               description of the required parameters associated with
%               an instrument.
%
%               [1,c1] = size(inst), char = class(inst)
%
%      raydir   the ray pointing toward a target. The ray emanates from
%               the location of the ephemeris object designated by the
%               input argument 'obsrvr' and is expressed relative to the
%               reference frame designated by 'rframe' (see descriptions
%               below).
%
%                [3,1] = size(raydir), double = class(raydir)
%
%      rframe   the string naming the reference frame associated with
%               the input ray's direction vector 'raydir'.
%
%               Since light time corrections are not supported for
%               rays, the orientation of the frame is always evaluated
%               at the epoch associated with the observer, as opposed
%               to the epoch associated with the light-time corrected
%               position of the frame center.
%
%               Case and leading or trailing blanks bracketing a non-blank
%               frame name are not significant in the string 'rframe'.
%
%               [1,c2] = size(rframe), char = class(rframe)
%
%      abcorr   the string indicating the aberration corrections
%               to apply when computing the 'raydir' direction.
%
%               The supported aberration correction options are
%
%                  'NONE'          No correction.
%                  'S'             Stellar aberration correction,
%                                  reception case.
%                  'XS'            Stellar aberration correction,
%                                  transmission case.
%
%               For detailed information, see the geometry finder required
%               reading, gf.req.
%
%               Case, leading and trailing blanks are not significant
%               in the string 'abcorr'.
%
%               [1,c3] = size(abcorr), char = class(abcorr)
%
%      obsrvr   the string naming the body from which the target
%               represented by 'raydir' is observed. The instrument designated
%               by 'inst' is treated as if it were co-located with the observer.
%
%               Optionally, you may supply the ID code of the object as an
%               integer string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the observer
%               is Earth.
%
%               [1,c4] = size(obsrvr), char = class(obsrvr)
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
%      nintvls  the maximum number of intervals to return in 'result'.
%               Note: this value should equal at least the number of expected
%               intervals. Recall two double precision values define
%               an interval.
%
%               [1,1] = size(nintvls), int32 = class(nintvls)
%
%   the call:
%
%      result = cspice_gfrfov( inst, raydir, rframe, abcorr, ...
%                              obsrvr, step, cnfine, nintvls)
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
%               [2n,1] = size(result), double = class(result)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      This example is an extension of the example in the
%      header of
%
%         cspice_gftfov
%
%      The problem statement for that example is
%
%
%         Search for times when Saturn's satellite Phoebe is within the
%         FOV of the Cassini narrow angle camera (CASSINI_ISS_NAC). To
%         simplify the problem, restrict the search to a short time
%         period where continuous Cassini bus attitude data are
%         available.
%
%         Use a step size of 10 seconds to reduce chances of missing
%         short visibility events.
%
%
%      Here we search the same confinement window for times when a
%      selected background star is visible. We use the FOV of the
%      Cassini ISS wide angle camera (CASSINI_ISS_WAC) to enhance the
%      probability of viewing the star.
%
%      The star we'll use has catalog number 6000 in the Hipparcos
%      Catalog. The star's J2000 right ascension and declination,
%      proper motion, and parallax are taken from that catalog.
%
%      Use the meta-kernel from the cspice_gftfov example:
%
%
%           KPL/MK
%
%           File name: gftfov_ex1.tm
%
%           This meta-kernel is intended to support operation of SPICE
%           example programs. The kernels shown here should not be
%           assumed to contain adequate or correct versions of data
%           required by SPICE-based user applications.
%
%           In order for an application to use this meta-kernel, the
%           kernels referenced here must be present in the user's
%           current working directory.
%
%           The names and contents of the kernels referenced
%           by this meta-kernel are as follows:
%
%              File name                     Contents
%              ---------                     --------
%              naif0009.tls                  Leapseconds
%              cpck05Mar2004.tpc             Satellite orientation and
%                                            radii
%              981005_PLTEPH-DE405S.bsp      Planetary ephemeris
%              020514_SE_SAT105.bsp          Satellite ephemeris
%              030201AP_SK_SM546_T45.bsp     Spacecraft ephemeris
%              cas_v37.tf                    Cassini FK
%              04135_04171pc_psiv2.bc        Cassini bus CK
%              cas00084.tsc                  Cassini SCLK kernel
%              cas_iss_v09.ti                Cassini IK
%
%
%           \begindata
%
%              KERNELS_TO_LOAD = ( 'naif0009.tls',
%                                  'cpck05Mar2004.tpc',
%                                  '981005_PLTEPH-DE405S.bsp',
%                                  '020514_SE_SAT105.bsp',
%                                  '030201AP_SK_SM546_T45.bsp',
%                                  'cas_v37.tf',
%                                  '04135_04171pc_psiv2.bc',
%                                  'cas00084.tsc',
%                                  'cas_iss_v09.ti' )
%           \begintext
%
%
%      MAXWIN  =  1000;
%      TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%      AU      =  149597870.693;
%
%      %
%      % Load kernels.
%      %
%      cspice_furnsh( 'gftfov_ex1.tm' )
%
%      %
%      % Store the time bounds of our search interval in
%      % the 'cnfine' confinement window.
%      %
%      et = cspice_str2et( { '2004 JUN 11 06:30:00 TDB', ...
%                            '2004 JUN 11 12:00:00 TDB' } );
%
%      cnfine = cspice_wninsd( et(1), et(2) );
%
%      %
%      %Initialize inputs for the search.
%      %
%      inst = 'CASSINI_ISS_WAC';
%
%      %
%      % Create a unit direction vector pointing from observer to star.
%      % We'll assume the direction is constant during the confinement
%      % window, and we'll use 'et[0]' as the epoch at which to compute the
%      % direction from the spacecraft to the star.
%      %
%      % The data below are for the star with catalog number 6000
%      % in the Hipparcos catalog. Angular units are degrees; epochs
%      % have units of Julian years and have a reference epoch of J1950.
%      % The reference frame is J2000.
%      %
%
%      parallax_deg = 0.000001056;
%
%      ra_deg_0     = 19.290789927;
%      ra_pm        = -0.000000720;
%      ra_epoch     = 41.2000;
%
%      dec_deg_0    =  2.015271007;
%      dec_pm       =  0.000001814;
%      dec_epoch    = 41.1300;
%
%      rframe       = 'J2000';
%      nintvls      = MAXWIN;
%
%
%      %
%      % Correct the star's direction for proper motion.
%      %
%      % The argument 't' represents 'et[0]' as Julian years past J1950.
%      %
%
%      t       = et(1)/cspice_jyear + ( cspice_j2000- cspice_j1950 )/365.25;
%
%      dtra    = t - ra_epoch;
%      dtdec   = t - dec_epoch;
%
%      ra_deg  = ra_deg_0  + dtra * ra_pm;
%      dec_deg = dec_deg_0 + dtra * dec_pm;
%
%      ra      = ra_deg  * cspice_rpd;
%      dec     = dec_deg * cspice_rpd;
%
%      starpos = cspice_radrec( 1, ra, dec );
%
%
%      %
%      % Correct star position for parallax applicable at
%      % the Cassini orbiter's position. (The parallax effect
%      % is negligible in this case; we're simply demonstrating
%      % the computation.)
%      %
%      parallax = parallax_deg * cspice_rpd;
%
%      stardist = AU/tan(parallax);
%
%      %
%      % Scale the star's direction vector by its distance from
%      % the solar system barycenter. Subtract off the position
%      % of the spacecraft relative to the solar system barycenter;
%      % the result is the ray's direction vector.
%      %
%      starpos = stardist * starpos;
%
%      [pos, ltime] = cspice_spkpos( 'cassini', et(1), 'J2000',  ...
%                                   'NONE', 'solar system barycenter' );
%
%      raydir = starpos - pos;
%
%      %
%      % Correct the star direction for stellar aberration when
%      % we conduct the search.
%      %
%      abcorr = 'S';
%      obsrvr = 'CASSINI';
%      stepsz = 100;
%
%      result = cspice_gfrfov( inst,   raydir, rframe, abcorr, ...
%                              obsrvr, stepsz, cnfine, nintvls );
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
%            disp( ' ' )
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
%      From : 2004-JUN-11 06:30:00.000000 (TDB)
%      To   : 2004-JUN-11 12:00:00.000000 (TDB)
%
%-Particulars
%
%   This routine determines a set of one or more time intervals when
%   the specified ray in contained within the field of view of a
%   specified instrument. We'll use the term "visibility event" to
%   designate such an appearance. The set of time intervals resulting
%   from the search is returned as a SPICE window.
%
%   The Search Process
%   ==================
%
%   The search for visibility events is treated as a search for state
%   transitions: times are sought when the state of the ray
%   changes from "not visible" to "visible" or vice versa.
%
%   Step Size
%   =========
%
%   Each interval of the confinement window is searched as follows:
%   first, the input step size is used to determine the time
%   separation at which the visibility state will be sampled.
%   Starting at the left endpoint of an interval, samples will be
%   taken at each step. If a state change is detected, a root has
%   been bracketed; at that point, the "root"--the time at which the
%   state change occurs---is found by a refinement process, for
%   example, by a binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the visibility state is constant:
%   the step size should be shorter than the shortest visibility event
%   duration and the shortest period between visibility events, within
%   the confinement window.
%
%   Having some knowledge of the relative geometry of the ray and
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
%   the CSPICE routine gfrfov_c.
%
%   MICE.REQ
%   CK.REQ
%   FRAMES.REQ
%   GF.REQ
%   KERNEL.REQ
%   NAIF_IDS.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%   WINDOWS.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 05-NOV-2013, EDW (JPL)
%
%      Corrected minor typos in header.
%
%      Renamed the argument 'size' to 'nintvls'. "size" is a Matlab function
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
%   GF ray in instrument FOV search
%
%-&

function [result] = cspice_gfrfov( inst,   raydir, rframe, abcorr, ...
                                  obsrvr, step,   cnfine, nintvls)

  switch nargin

     case 8

        inst    = zzmice_str(inst);
        raydir  = zzmice_dp(raydir);
        rframe  = zzmice_str(rframe);
        abcorr  = zzmice_str(abcorr);
        obsrvr  = zzmice_str(obsrvr);
        step    = zzmice_dp(step);
        cnfine  = zzmice_win(cnfine);
        nintvls = zzmice_int(nintvls,    [1, int32(inf)/2] );

     otherwise

        error ( [ 'Usage: [result] = cspice_gfrfov( `inst`, ' ...
                               'raydir[3], `rframe` '         ...
                               '`abcorr`, `obsrvr`, step, cnfine, nintvls)' ] )

  end

  %
  % Call the GF routine, add to 'cnfine' the space needed for
  % the control segment.
  %
  try

     [result] = mice('gfrfov_c',  inst,   ...
                                  raydir, ...
                                  rframe, ...
                                  abcorr, ...
                                  obsrvr, ...
                                  step,   ...
                                  [zeros(6,1); cnfine], nintvls);

  catch
     rethrow(lasterror)
  end
