%-Abstract
%
%   CSPICE_GFILUM determines the time intervals over which a specified
%   constraint on the observed phase, solar incidence, or emission angle
%   at a specified target body surface point is met.
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
%      method    name specifying the computation method to use. Parameters
%                include, but are not limited to, the shape model
%                used to represent the surface of the target body.
%
%                The only choice currently supported is
%
%                   'Ellipsoid'        The illumination angle
%                                      computation uses a triaxial
%                                      ellipsoid to model the surface
%                                      of the target body. The
%                                      ellipsoid's radii must be
%                                      available in the kernel pool.
%
%                Neither case nor white space are significant in
%                'method'. For example, the string ' eLLipsoid ' is
%                valid.
%
%                [1,c1] = size(method); char = class(method)
%
%                  or
%
%                [1,1] = size(method); cell = class(method)
%
%      angtyp    name specifying the illumination angle for which a search
%                is to be performed. The possible values of 'angtyp' are:
%
%                   'PHASE'
%                   'SOLAR INCIDENCE'
%                   'EMISSION'
%
%                See the Particulars section below for a detailed
%                description of these angles.
%
%                Neither case nor white space are significant in
%                'angtyp'. For example, the string ' Solar incidence ' is
%                valid.
%
%                [1,c2] = size(angtyp); char = class(angtyp)
%
%                  or
%
%                [1,1] = size(angtyp); cell = class(angtyp)
%
%      target    name of the target body. The point at which the
%                illumination angles are defined is located on the
%                surface of this body.
%
%                Optionally, you may supply the integer ID code for
%                the object as an integer string. For example both
%                'MOON' and '301' are legitimate strings that indicate
%                the moon is the target body.
%
%                [1,c3] = size(target); char = class(target)
%
%                  or
%
%                [1,1] = size(target); cell = class(target)
%
%      illmn     name of the illumination source. This source
%                may be any ephemeris object. Case, blanks, and
%                numeric values are treated in the same way as for the
%                input 'target'.
%
%                [1,c4] = size(illmn); char = class(illmn)
%
%                  or
%
%                [1,1] = size(illmn); cell = class(illmn)
%
%      fixref    name of the body-fixed, body-centered
%                reference frame associated with the target body. The
%                input surface point 'spoint' is expressed relative to
%                this reference frame, and this frame is used to
%                define the orientation of the target body as a
%                function of time.
%
%                The string 'fixref' is case-insensitive, and leading
%                and trailing blanks in 'fixref' are not significant.
%
%                [1,c5] = size(fixref); char = class(fixref)
%
%                  or
%
%                [1,1] = size(fixref); cell = class(fixref)
%
%      abcorr    describes the aberration corrections to apply to the state
%                evaluations to account for one-way light time and stellar
%                aberration.
%
%                Any 'reception' correction accepted by spkezr_c can be
%                used here. See the header of spkezr_c for a detailed
%                description of the aberration correction options. For
%                convenience, the options are listed below:
%
%                   'NONE'     Apply no correction.
%
%                   'LT'       'Reception' case:  correct for
%                              one-way light time using a Newtonian
%                              formulation.
%
%                   'LT+S'     'Reception' case:  correct for
%                              one-way light time and stellar
%                              aberration using a Newtonian
%                              formulation.
%
%                   'CN'       'Reception' case:  converged
%                              Newtonian light time correction.
%
%                   'CN+S'     'Reception' case:  converged
%                              Newtonian light time and stellar
%                              aberration corrections.
%
%                Case and blanks are not significant in the string
%                'abcorr'.
%
%                [1,c6] = size(abcorr); char = class(abcorr)
%
%                  or
%
%                [1,1] = size(abcorr); cell = class(abcorr)
%
%      obsrvr    name of the observing body. Optionally, you
%                may supply the ID code of the object as an integer
%                string. For example, both 'EARTH' and '399' are
%                legitimate strings to supply to indicate that the
%                observer is Earth.
%
%                [1,c7] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%                [1,1] = size(obsrvr); cell = class(obsrvr)
%
%      spoint    a surface point on the target body, expressed in
%                Cartesian coordinates, relative to the body-fixed
%                target frame designated by 'fixref'.
%
%                'spoint' need not be visible from the observer's
%                location in order for the constraint specified by
%                'relate' and 'refval' (see descriptions below) to be
%                satisfied.
%
%                The components of 'spoint' have units of km.
%
%                [3,1] = size(spoint); double = class(spoint)
%
%      relate    describes the constraint relational operator on a specified
%                illumination angle. The result window found by this routine
%                indicates the time intervals where the constraint is
%                satisfied. Supported values of 'relate' and corresponding
%                meanings are shown below:
%
%                    '>'      The angle is greater than the reference
%                           value 'refval'.
%
%                    '='      The angle is equal to the reference
%                             value 'refval'.
%
%                    '<'      The angle is less than the reference
%                             value 'refval'.
%
%                   'ABSMAX'  The angle is at an absolute maximum.
%
%                   'ABSMIN'  The angle is at an absolute  minimum.
%
%                   'LOCMAX'  The angle is at a local maximum.
%
%                   'LOCMIN'  The angle is at a local minimum.
%
%                The caller may indicate that the region of interest is
%                the set of time intervals where the angle is within a
%                specified separation from an absolute extremum. The
%                argument 'adjust' (described below) is used to specify
%                this separation.
%
%                Local extrema are considered to exist only in the
%                interiors of the intervals comprising the confinement
%                window:  a local extremum cannot exist at a boundary
%                point of the confinement window.
%
%                Case is not significant in the string 'relate'.
%
%                [1,c8] = size(relate); char = class(relate)
%
%                  or
%
%                [1,1] = size(relate); cell = class(relate)
%
%      refval    reference value used together with the argument
%                'relate' to define an equality or inequality to be
%                satisfied by the specified illumination angle. See the
%                discussion of 'relate' above for further information.
%
%                The units of 'refval' are radians.
%
%                [1,1] = size(refval); double = class(refval)
%
%      adjust    parameter used to modify searches for absolute
%                extrema: when 'relate' is set to 'ABSMAX' or 'ABSMIN'
%                and 'adjust' is set to a positive value, gfilum_c will
%                find times when the observer-target distance is within
%                'adjust' km of the specified extreme value.
%
%                If 'adjust' is non-zero and a search for an absolute
%                minimum 'min' is performed, the result window contains
%                time intervals when the observer-target distance has
%                values between 'min' and min+adjust.
%
%                If the search is for an absolute maximum 'max', the
%                corresponding range is from max-adjust to 'max'.
%
%                'adjust' is not used for searches for local extrema,
%                equality or inequality conditions.
%
%               [1,1] = size(adjust); double = class(adjust)
%
%      step      step size to use in the search. 'step' must
%                be short enough for a search using this step size to
%                locate the time intervals where the specified
%                illumination angle is monotone increasing or
%                decreasing. However, 'step' must not be *too* short, or
%                the search will take an unreasonable amount of time.
%
%                The choice of 'step' affects the completeness but not
%                the precision of solutions found by this routine; the
%                precision is controlled by the convergence tolerance.
%                See the discussion of the parameter SPICE_GF_CNVTOL for
%                details.
%
%                'step' has units of seconds.
%
%                [1,1] = size(step); double = class(step)
%
%      nintvls   is a parameter specifying the number of intervals that
%                can be accommodated by each of the dynamically allocated
%                workspace windows used internally by this routine.
%
%                In many cases, it's not necessary to compute an accurate
%                estimate of how many intervals are needed; rather, the
%                user can pick a size considerably larger than what's
%                really required.
%
%                However, since excessively large arrays can prevent
%                applications from compiling, linking, or running
%                properly, sometimes 'nintvls' must be set according to
%                the actual workspace requirement. A rule of thumb for
%                the number of intervals needed is
%
%                   nintvls  =  2*n  +  ( m / step )
%
%                where
%
%                   n     is the number of intervals in the confinement
%                         window
%
%                   m     is the measure of the confinement window, in
%                         units of seconds
%
%                   step  is the search step size in seconds
%
%                [1,1] = size(nintvls); int32 = class(nintvls)
%
%      cnfine    a SPICE window that confines the time period over
%                which the specified search is conducted. 'cnfine' may
%                consist of a single interval or a collection of
%                intervals.
%
%                The endpoints of the time intervals comprising 'cnfine'
%                are interpreted as seconds past J2000 TDB.
%
%                See the Examples section below for a code example that
%                shows how to create a confinement window.
%
%                [2r,1] = size(cnfine), double = class(cnfine)
%
%   the call:
%
%      result = cspice_gfilum( method, angtyp,  target, illmn,  ...
%                              fixref, abcorr,  obsrvr, spoint, ...
%                              relate, refval,  adjust,        ...
%                              step,   nintvls, cnfine )
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
%               [2s,1] = size(result), double = class(result)
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
%
%         KPL/MK
%
%         File name: gfilum.tm
%
%         This is the meta-kernel file for the example problem for
%         the subroutine gfilum_c. These kernel files can be found on
%         the NAIF website.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                             Contents
%            ---------                             --------
%            de421.bsp                             Planetary ephemeris
%            pck00010.tpc                          Planet orientation
%                                                  and radii
%            naif0010.tls                          Leapseconds
%            spk_psp_110101-130101_081216_3pm.bsp  MRO predict SPK
%
%
%         \begindata
%
%         KERNELS_TO_LOAD = ( 'de421.bsp',
%                             'pck00010.tpc',
%                             'naif0010.tls',
%                             'spk_psp_110101-130101+'
%                             '_081216_3pm.bsp'       )
%         \begintext
%
%         End of meta-kernel
%
%   Example:
%
%      Determine time intervals over which the planned Mars Science
%      Laboratory (MSL) Gale Crater landing site satisfies certain
%      constraints on its illumination and visibility as seen from
%      the Mars Reconnaissance Orbiter (MRO) spacecraft. The
%      observation period will range from slightly before the planned
%      landing time to about 10 days later.
%
%      In this case we require the emission angle to be less than
%      30 degrees and the solar incidence angle to be less than
%      40 degrees.
%
%      %
%      % Output time format:
%      %
%      TIMFMT = 'YYYY MON DD HR:MN:SC.###### TDB::TDB';
%
%      %
%      % Meta-kernel name:
%      %
%      META = 'gfilum.tm';
%
%      %
%      % Maximum number of intervals in the windows
%      % used in this program:
%      %
%      MAXIVL = 1000;
%
%      %
%      % Local variables
%      %
%      r2d    = cspice_dpr();
%
%      %
%      % Initial values
%      %
%      % Mars planetodetic coordinates of landing site.
%      % Angular units are degrees; distance units are km.
%      %
%      gclat  =  -4.543182;
%      gclon  = 137.420000;
%      gcalt  =  -4.876405;
%
%      %
%      % Load kernels:
%      %
%      cspice_furnsh( META );
%
%      %
%      % Convert the landing site location from planetodetic
%      % to Cartesian coordinates for use with GFILUM.
%      %
%      radii = cspice_bodvrd( 'MARS', 'RADII', 3 );
%
%      re = radii(1);
%      rp = radii(3);
%
%      f  = ( re - rp ) / re;
%
%      gcpos = cspice_georec( gclon * cspice_rpd(), ...
%                             gclat * cspice_rpd(), ...
%                             gcalt, re, f);
%
%      %
%      % Set the search interval:
%      %
%      utcbeg = '2012 AUG 5 00:00:00 UTC';
%      et0    = cspice_str2et( utcbeg );
%
%      utcend = '2012 SEP 15 00:00:00 UTC';
%      et1    = cspice_str2et( utcend );
%
%      cnfine = cspice_wninsd( et0, et1 );
%
%
%      %
%      % Set observer, target, aberration correction, and the
%      % Mars body-fixed, body-centered reference frame. The
%      % lighting source is the sun.
%      %
%      % Aberration corrections are set for remote observations.
%      %
%      illmn  = 'sun';
%      obsrvr = 'mro';
%      target = 'mars';
%      abcorr = 'cn+s';
%      fixref = 'iau_mars';
%
%      %
%      % Initialize the adjustment value for absolute
%      % extremum searches. We're not performing
%      % such searches in this example, but this input
%      % to GFILUM must still be set.
%      %
%      adjust = 0.0;
%
%      %
%      % The computation uses an ellipsoidal model for the
%      % target body shape.
%      %
%      method = 'Ellipsoid';
%
%      %
%      % Set the reference value to use for the solar
%      % incidence angle search.
%      %
%      refval = 45.0 * cspice_rpd();
%
%      %
%      % Since the period of the solar incidence angle
%      % is about one Martian day, we can safely use 6 hours
%      % as the search step.
%      %
%      step   = 21600.0;
%
%      %
%      % Search over the confinement window for times
%      % when the solar incidence angle is less than
%      % the reference value.
%      %
%      [wnsolr] = cspice_gfilum( method, 'SOLAR INCIDENCE', target, ...
%                                illmn,  fixref,            abcorr, ...
%                                obsrvr, gcpos,             '<',    ...
%                                refval, adjust,            step,   ...
%                                MAXIVL, cnfine );
%
%      %
%      % With the search on the incidence angle complete, perform
%      % a search on the emission angle.
%      %
%      % Set the reference value for the emission angle search.
%      %
%      refval = 80.0 * cspice_rpd();
%
%      %
%      % We'll use 15 minutes as the search step. This step
%      % is small enough to be suitable for Mars orbiters.
%      % Units are seconds.
%      %
%      step   = 900.0;
%
%      %
%      % Search over the previous result window for times when the
%      % emission angle is less than the reference value.
%      %
%      [result] = cspice_gfilum( method, 'EMISSION', target, illmn, ...
%                                fixref, abcorr,     obsrvr, gcpos, ...
%                                '<',    refval,     adjust, step,  ...
%                                MAXIVL, wnsolr );
%
%      %
%      % Display the result window. Show the solar incidence
%      % and emission angles at the window's interval
%      % boundaries.
%      %
%      if ( cspice_wncard( result ) == 0 )
%
%         disp( '     Window is empty: condition is not met.' )
%
%      else
%
%         fprintf( '                                     ' )
%         fprintf( '       Solar Incidence   Emission\n'   )
%         fprintf( '                                     ' )
%         fprintf( '             (deg)         (deg)\n\n'  )
%
%         for i=1:cspice_wncard( result )
%
%            [start, finish] = cspice_wnfetd( result, i );
%
%            %
%            % Compute the angles of interest at the boundary
%            % epochs.
%            %
%            timstr = cspice_timout( start, TIMFMT );
%            [trgepc, srfvec, phase, solar, emissn] =                ...
%                                     cspice_ilumin( method, target, ...
%                                                    start,  fixref, ...
%                                                    abcorr, obsrvr, ...
%                                                    gcpos );
%
%               fprintf ( '    Start: %s %14.9f %14.9f\n', timstr, ...
%                                                       solar*r2d, ...
%                                                       emissn*r2d )
%
%
%            timstr = cspice_timout( finish, TIMFMT);
%            [trgepc, srfvec, phase, solar, emissn] =                ...
%                                     cspice_ilumin( method, target, ...
%                                                    finish, fixref, ...
%                                                    abcorr, obsrvr, ...
%                                                    gcpos );
%
%               fprintf ( '    Start: %s %14.9f %14.9f\n\n', timstr,...
%                                                       solar*r2d,  ...
%                                                       emissn*r2d )
%
%         end
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%                                              Solar Incidence   Emission
%                                                    (deg)         (deg)
%
%      Start: 2012 AUG 12 08:20:03.526171 TDB   43.162554493   80.000000005
%      Start: 2012 AUG 12 08:24:01.612468 TDB   44.061827678   80.000000003
%
%      Start: 2012 AUG 17 11:45:09.843860 TDB   44.655209578   79.999999988
%      Start: 2012 AUG 17 11:46:39.910178 TDB   45.000000000   75.190270778
%
%      Start: 2012 AUG 23 15:32:01.989001 TDB   42.039781850   80.000000002
%      Start: 2012 AUG 23 15:33:58.603546 TDB   42.488410630   80.000000002
%
%      Start: 2012 AUG 28 18:56:36.440913 TDB   43.450078325   80.000000007
%      Start: 2012 AUG 28 19:01:25.343034 TDB   44.575503232   79.999999990
%
%      Start: 2012 SEP 09 02:08:12.664420 TDB   42.294339879   80.000000007
%      Start: 2012 SEP 09 02:11:45.039261 TDB   43.133576379   80.000000004
%
%      Start: 2012 SEP 14 05:33:12.767223 TDB   43.878824441   80.000000018
%      Start: 2012 SEP 14 05:37:54.330789 TDB   45.000000001   77.289478337
%
%-Particulars
%
%   This routine determines a set of one or more time intervals
%   within the confinement window when the specified illumination
%   angle satisfies a caller-specified constraint. The resulting set
%   of intervals is returned as a SPICE window.
%
%   The term 'illumination angles' refers to the following set of
%   angles:
%
%
%      phase angle              Angle between the vectors from the
%                               surface point to the observer and
%                               from the surface point to the
%                               illumination source.
%
%      incidence angle          Angle between the surface normal at
%                               the specified surface point and the
%                               vector from the surface point to the
%                               illumination source.
%
%      emission angle           Angle between the surface normal at
%                               the specified surface point and the
%                               vector from the surface point to the
%                               observer.
%
%   The diagram below illustrates the geometric relationships
%   defining these angles. The labels for the incidence, emission,
%   and phase angles are 'inc.', 'e.', and 'phase'.
%
%
%                                                    *
%                                            illumination source
%
%                  surface normal vector
%                            ._                 _.
%                            |\                 /|  illumination
%                              \    phase      /    source vector
%                               \   .    .    /
%                               .            .
%                                 \   ___   /
%                            .     \/     \/
%                                  _\ inc./
%                           .    /   \   /
%                           .   |  e. \ /
%       *             <--------------- *  surface point on
%    viewing            vector            target body
%    location           to viewing
%    (observer)         location
%
%
%   Note that if the target-observer vector, the target normal vector
%   at the surface point, and the target-illumination source vector
%   are coplanar, then phase is the sum of the incidence and emission
%   angles. This rarely occurs; usually
%
%      phase angle  <  incidence angle + emission angle
%
%   All of the above angles can be computed using light time
%   corrections, light time and stellar aberration corrections, or no
%   aberration corrections. In order to describe apparent geometry as
%   observed by a remote sensing instrument, both light time and
%   stellar aberration corrections should be used.
%
%   The way aberration corrections are applied by this routine
%   is described below.
%
%      Light time corrections
%      ======================
%
%         Observer-target surface point vector
%         ------------------------------------
%
%         Let ET be the epoch at which an observation or remote
%         sensing measurement is made, and let ET - LT ('LT' stands
%         for 'light time') be the epoch at which the photons
%         received at ET were emitted from the surface point 'spoint'.
%         Note that the light time between the surface point and
%         observer will generally differ from the light time between
%         the target body's center and the observer.
%
%
%         Target body's orientation
%         -------------------------
%
%         Using the definitions of ET and LT above, the target body's
%         orientation at ET - LT is used. The surface normal is
%         dependent on the target body's orientation, so the body's
%         orientation model must be evaluated for the correct epoch.
%
%
%         Target body -- illumination source vector
%         -----------------------------------------
%
%         The surface features on the target body near 'spoint' will
%         appear in a measurement made at ET as they were at ET-LT.
%         In particular, lighting on the target body is dependent on
%         the apparent location of the illumination source as seen
%         from the target body at ET-LT. So, a second light time
%         correction is used to compute the position of the
%         illumination source relative to the surface point.
%
%
%      Stellar aberration corrections
%      ==============================
%
%      Stellar aberration corrections are applied only if
%      light time corrections are applied as well.
%
%         Observer-target surface point body vector
%         -----------------------------------------
%
%         When stellar aberration correction is performed, the
%         observer-to-surface point direction vector, which we'll
%         call SRFVEC, is adjusted so as to point to the apparent
%         position of 'spoint': considering 'spoint' to be an ephemeris
%         object, SRFVEC points from the observer's position at ET to
%         the light time and stellar aberration
%         corrected position of 'spoint'.
%
%         Target body-illumination source vector
%         --------------------------------------
%
%         The target body-illumination source vector is the apparent
%         position of the illumination source, corrected for light
%         time and stellar aberration, as seen from the surface point
%         'spoint' at time ET-LT.
%
%
%   Below we discuss in greater detail aspects of this routine's
%   solution process that are relevant to correct and efficient
%   use of this routine in user applications.
%
%
%   The Search Process
%   ==================
%
%   Regardless of the type of constraint selected by the caller, this
%   routine starts the search for solutions by determining the time
%   periods, within the confinement window, over which the specified
%   illumination angle is monotone increasing and monotone decreasing.
%   Each of these time periods is represented by a SPICE window.
%   Having found these windows, all of the illumination angle's local
%   extrema within the confinement window are known. Absolute extrema
%   then can be found very easily.
%
%   Within any interval of these 'monotone' windows, there will be at
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
%
%   Step Size
%   =========
%
%   The monotone windows (described above) are found via a two-step
%   search process. Each interval of the confinement window is
%   searched as follows: first, the input step size is used to
%   determine the time separation at which the sign of the rate of
%   change of the illumination angle will be sampled. Starting at
%   the left endpoint of an interval, samples will be taken at each
%   step. If a change of sign is found, a root has been bracketed; at
%   that point, the time at which the range rate is zero can be
%   found by a refinement process, for example, via binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the illumination angle is monotone:
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
%
%   Convergence Tolerance
%   =====================
%
%   As described above, the root-finding process used by this routine
%   involves first bracketing roots and then using a search process
%   to locate them. 'Roots' are both times when local extrema are
%   attained and times when the illumination angle is equal to a
%   reference value. All endpoints of the intervals comprising the
%   result window are either endpoints of intervals of the
%   confinement window or roots.
%
%   Once a root has been bracketed, a refinement process is used to
%   narrow down the time interval within which the root must lie.
%   This refinement process terminates when the location of the root
%   has been determined to within an error margin called the
%   'convergence tolerance.' The convergence tolerance used by this
%   routine is set via the parameter SPICE_GF_CNVTOL.
%
%   The value of SPICE_GF_CNVTOL is set to a 'tight' value so that the
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
%   the CSPICE routine gfilum_c.
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
%   -Mice Version 1.0.0, 07-NOV-2013, EDW (JPL)
%
%-Index_Entries
%
%   solve for illumination_angle constraints
%   solve for phase_angle constraints
%   solve for solar_incidence_angle constraints
%   solve for emission_angle constraints
%   search using illumination_angle constraints
%
%-&

function [result] = cspice_gfilum( method, angtyp,  target, illmn,  ...
                                   fixref, abcorr,  obsrvr, spoint, ...
                                   relate, refval,  adjust,         ...
                                   step,   nintvls, cnfine )

   switch nargin

      case 14

         method  = zzmice_str(method);
         angtyp  = zzmice_str(angtyp);
         target  = zzmice_str(target);
         illmn   = zzmice_str(illmn);
         fixref  = zzmice_str(fixref);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         spoint  = zzmice_dp(spoint);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);

      otherwise

         error ( [ 'Usage: [result] = cspice_gfilum( `method`, `angtyp`, '...
                                    '`target`, `illmn`, `fixref`, '       ...
                                    '`abcorr`, `obsrvr`, spoint[3], '     ...
                                    '`relate`, refval, adjust, '          ...
                                    'step, nintvls, cnfine, )' ] )

   end

   %
   % Call the GF routine, add to 'cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gfilum_c', method,  angtyp, target, illmn,  ...
                                  fixref,  abcorr, obsrvr, spoint, ...
                                  relate,  refval, adjust, step,   ...
                                  nintvls, [zeros(6,1); cnfine] );
   catch
      rethrow(lasterror)
   end
