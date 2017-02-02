%-Abstract
%
%   CSPICE_GFTFOV determines time intervals when a specified ephemeris
%   object intersects the space bounded by the field-of-view (FOV) of a
%   specified instrument.
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
%               (FOV) of which is to be used for a target intersection
%               search: times when the specified target intersects the
%               region of space corresponding to the FOV are sought.
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
%      target   the string naming the 'target' body, the appearances
%               of which in the specified instrument's field of view are
%               sought. The body must be an ephemeris object.
%
%               Optionally, you may supply the integer NAIF ID code
%               for the body as a string. For example both 'MOON' and
%               '301' are legitimate strings that designate the Moon.
%
%               The 'target' string lacks sensitivity to case, and to leading
%               and trailing blanks.
%
%               [1,c2] = size(target), char = class(target)
%
%      tshape   the string naming the geometric model used to
%               represent the shape of the 'target' body. The supported
%               options are:
%
%                  'ELLIPSOID'   Use a triaxial ellipsoid model,
%                                with radius values provided via the
%                                kernel pool. A kernel variable
%                                having a name of the form
%
%                                   'BODYnnn_RADII'
%
%                                where nnn represents the NAIF
%                                integer code associated with the
%                                body, must be present in the kernel
%                                pool. This variable must be
%                                associated with three numeric
%                                values giving the lengths of the
%                                ellipsoid's X, Y, and Z semi-axes.
%
%                  'POINT'       Treat the body as a single point.
%
%               The 'tshape' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%               [1,c3] = size(tshape), char = class(tshape)
%
%      tframe   the string naming the body-fixed, body-centered
%               reference frame associated with the target body. Examples of
%               such names are 'IAU_SATURN' (for Saturn) and 'ITRF93'
%               (for the Earth).
%
%               If the target body is modeled as a point, 'tframe'
%               is ignored and should be left blank.
%
%               The 'tframe' string lacks sensitivity to case, and to leading
%               and trailing blanks.
%
%               [1,c4] = size(tframe), char = class(tframe)
%
%      abcorr   the string indicating the aberration corrections to apply
%               to the state evaluations to account for one-way light time and
%               stellar aberration.
%
%               For remote sensing applications, where the apparent
%               position and orientation of the target seen by the
%               observer are desired, normally either of the
%               corrections
%
%                  'LT+S'
%                  'CN+S'
%
%               should be used. These and the other supported options
%               are described below.
%
%                 'NONE'      Apply no correction.
%
%               Supported aberration correction options for reception case
%               (radiation is received by observer at ET) are:
%
%                  'LT'       Correct for one-way light time using a Newtonian
%                             formulation.
%
%                  'LT+S'     Correct for one-way light time and stellar
%                             aberration using a Newtonian formulation.
%
%                  'CN'       Correct for one-way light time using a converged
%                             Newtonian light time correction.
%
%                  'CN+S'     Correct for one-way light time and stellar
%                             aberration using a converged Newtonian light time
%                             and stellar aberration corrections.
%
%               Supported aberration correction options for transmission case
%               (radiation is emitted from observer at ET) are:
%
%                  'XLT'      Correct for one-way light time using a Newtonian
%                             formulation.
%
%                  'XLT+S'    Correct for one-way light time and stellar
%                             aberration using a Newtonian formulation.
%
%                  'XCN'      Correct for one-way light time using a converged
%                             Newtonian light time correction.
%
%                  'XCN+S'    Correct for one-way light time and stellar
%                             aberration using a converged Newtonian light time
%                             and stellar aberration corrections.
%
%               For detailed information, see the geometry finder
%               required reading, gf.req.
%
%               The 'abcorr' string lacks sensitivity to case, and to leading
%               and trailing blanks.
%
%               [1,c5] = size(abcorr), char = class(abcorr)
%
%      obsrvr   the string naming the body from which the target is
%               observed. The instrument designated by 'inst' is treated
%               as if it were co-located with the observer.
%
%               Optionally, you may supply the ID code of the object as an
%               integer string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the observer
%               is Earth.
%
%               [1,c6] = size(abcorr), char = class(abcorr)
%
%      step     the step size to use in the search to use in the search.
%               'step' must be short enough for a search using step
%               to locate the time intervals where the specified
%               angular separation function is monotone increasing or
%               decreasing. However, 'step' must not be *too* short, or
%               the search will take an unreasonable amount of time.
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
%      result = cspice_gftfov( inst,   target, tshape, tframe, ...
%                              abcorr, obsrvr, step,   cnfine, room)
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
%      Search for times when Saturn's satellite Phoebe is within
%      the FOV of the Cassini narrow angle camera (CASSINI_ISS_NAC).
%      To simplify the problem, restrict the search to a short time
%      period where continuous Cassini bus attitude data are
%      available.
%
%      Use a step size of 10 seconds to reduce chances of missing
%      short visibility events.
%
%      Use the meta-kernel listed below to load the required SPICE
%      kernels.
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
%
%      %
%      % Load kernels.
%      %
%      cspice_furnsh( 'gftfov_ex1.tm' )
%
%      %
%      % Store the time bounds of our search interval in
%      % the cnfine confinement window.
%      %
%      et = cspice_str2et( { '2004 JUN 11 06:30:00 TDB', ...
%                            '2004 JUN 11 12:00:00 TDB' } );
%
%      cnfine = cspice_wninsd( et(1), et(2) );
%
%      %
%      %Initialize inputs for the search.
%      %
%      inst   = 'CASSINI_ISS_NAC';
%      target = 'PHOEBE';
%      tshape = 'ELLIPSOID';
%      tframe = 'IAU_PHOEBE';
%      abcorr = 'LT+S';
%      obsrvr = 'CASSINI';
%      step   = 10.;
%      room   = MAXWIN;
%
%      result = cspice_gftfov( inst,  target, tshape, tframe, ...
%                              abcorr, obsrvr, step, cnfine, room);
%
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
%      From : 2004-JUN-11 07:35:49.958590 (TDB)
%      To   : 2004-JUN-11 08:48:27.485966 (TDB)
%
%      From : 2004-JUN-11 09:03:19.767800 (TDB)
%      To   : 2004-JUN-11 09:35:27.634791 (TDB)
%
%      From : 2004-JUN-11 09:50:19.585474 (TDB)
%      To   : 2004-JUN-11 10:22:27.854254 (TDB)
%
%      From : 2004-JUN-11 10:37:19.332697 (TDB)
%      To   : 2004-JUN-11 11:09:28.116017 (TDB)
%
%      From : 2004-JUN-11 11:24:19.049485 (TDB)
%      To   : 2004-JUN-11 11:56:28.380305 (TDB)
%
%-Particulars
%
%   This routine determines a set of one or more time intervals
%   within the confinement window when any portion of a specified
%   target body appears within the field of view of a specified
%   instrument. We'll use the term "visibility event" to designate
%   such an appearance. The set of time intervals resulting from the
%   search is returned as a SPICE window.
%
%   Below we discuss in greater detail aspects of this routine's
%   solution process that are relevant to correct and efficient use
%   of this routine in user applications.
%
%   The Search Process
%   ==================
%
%   The search for visibility events is treated as a search for state
%   transitions: times are sought when the state of the target body
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
%   example, via binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the visibility state is constant:
%   the step size should be shorter than the shortest visibility event
%   duration and the shortest period between visibility events, within
%   the confinement window.
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
%   the CSPICE routine gftfov_c.
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
%   -Mice Version 1.1.0, 12-MAY-2012, EDW (JPL)
%
%      Corrected minor typo in header.
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
%   GF target in instrument FOV search
%
%-&

function [result] = cspice_gftfov( inst,  target, tshape, tframe, ...
                                  abcorr, obsrvr, step, cnfine, room)

   switch nargin

      case 9

         inst    = zzmice_str(inst);
         target  = zzmice_str(target);
         tshape  = zzmice_str(tshape);
         tframe  = zzmice_str(tframe);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         step    = zzmice_dp(step);
         cnfine  = zzmice_win(cnfine);
         room    = zzmice_int(room,    [1, int32(inf)/2] );

      otherwise

         error ( [ 'Usage: [result] = cspice_gftfov( `inst`, ' ...
                                   '`target`, `tshape`, `tframe` ' ...
                                   '`abcorr`, `obsrvr`, step, '  ...
                                   'cnfine, room)' ] )

   end

   %
   % Call the GF routine, add to 'cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gftfov_c',  inst,   ...
                                   target, ...
                                   tshape, ...
                                   tframe, ...
                                   abcorr, ...
                                   obsrvr, ...
                                   step,   ...
                                   [zeros(6,1); cnfine], room);

   catch
      rethrow(lasterror)
   end





