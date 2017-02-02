%-Abstract
%
%   CSPICE_EDTERM computes a set of points on the umbral or penumbral
%   terminator of a specified target body, where the target shape is modeled
%   as an ellipsoid.
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
%      trmtyp   string indicating the type of terminator to
%               compute: umbral or penumbral. The umbral terminator is
%               the boundary of the portion of the ellipsoid surface in
%               total shadow. The penumbral terminator is the boundary
%               of the portion of the surface that is completely
%               illuminated. Note that in astronomy references, the
%               unqualified word "terminator" refers to the umbral
%               terminator. Here, the unqualified word refers to either
%               type of terminator.
%
%               Possible values of 'trmtyp' are
%
%                  'UMBRAL'
%                  'PENUMBRAL'
%
%               Case and leading or trailing blanks in 'trmtyp' are
%               not significant.
%
%               [1,c1] = size(trmtyp); char = class(trmtyp)
%
%      source   string name of the body acting as a light source.
%               'source' is case-insensitive, and leading and trailing
%               blanks in 'target' are not significant. Optionally, you
%               may supply a string containing the integer ID code for
%               the object. For example both "SUN" and "10" are
%               legitimate strings that indicate the Sun is the light
%               source.
%
%               This routine assumes that a kernel variable representing
%               the light source's radii is present in the kernel pool.
%               Normally the kernel variable would be defined by loading
%               a PCK file.
%
%               The shape of the light source is always modeled as a
%               sphere, regardless of whether radii defining a triaxial
%               ellipsoidal shape model are available in the kernel
%               pool. The maximum radius of the body is used as the
%               radius of the sphere.
%
%               [1,c2] = size(source); char = class(source)
%
%      target   string name of the target body. 'target' is
%               case-insensitive, and leading and trailing blanks in
%               'target' are not significant. Optionally, you may supply
%               a string containing the integer ID code for the object.
%               For example both "MOON" and "301" are legitimate strings
%               that indicate the moon is the target body.
%
%               This routine assumes that a kernel variable representing
%               the target's radii is present in the kernel pool.
%               Normally the kernel variable would be defined by loading
%               a PCK file.
%
%               [1,c3] = size(target); char = class(target)
%
%      et       epoch of participation of the observer, expressed
%               as ephemeris seconds past J2000 TDB: 'et' is the epoch
%               at which the observer's position is computed.
%
%               When aberration corrections are not used, 'et' is also
%               the epoch at which the position and orientation of the
%               target body and position of the light source are
%               computed.
%
%               When aberration corrections are used, 'et' is the epoch
%               at which the observer's position relative to the solar
%               system barycenter is computed; in this case the position
%               and orientation of the target body are computed at
%               et-lt, where lt is the one-way light time between the
%               target body's center and the observer. See the
%               description of 'abcorr' below for details.
%
%               [1,1] = size(et); double = class(et)
%
%      fixref   string name of the reference frame relative to
%               which the output terminator points are expressed. This must
%               be a body-centered, body-fixed frame associated with the
%               target. The frame's axes must be compatible with the
%               triaxial ellipsoidal shape model associated with the
%               target body (normally provide via a PCK): this routine
%               assumes that the first, second, and third axis lengths
%               correspond, respectively, to the x, y, and z-axes of the
%               frame designated by 'fixref'.
%
%               'fixref' may refer to a built-in frame (documented in
%               the Frames Required Reading) or a frame defined by a
%               loaded frame kernel (FK).
%
%               The orientation of the frame designated by 'fixref' is
%               evaluated at epoch of participation of the target body.
%               See the descriptions of 'et' and 'abcorr' for details.
%
%               [1,c4] = size(fixref); char = class(fixref)
%
%      abcorr   string indicating the aberration correction to be
%               applied when computing the observer-target position, the
%               orientation of the target body, and the target-
%               source position vector. 'abcorr' may be any of
%               the following.
%
%                  'NONE'     Apply no correction. Compute the
%                             terminator points using the position
%                             of the light source and target, and
%                             the orientation of the target, at 'et'.
%
%               Let 'lt' represent the one-way light time between the
%               observer and the target body's center. The following
%               values of 'abcorr' apply to the "reception" case in
%               which photons depart from the target body's center at
%               the light-time corrected epoch et-lt and *arrive* at
%               the observer's location at 'et':
%
%
%                  'LT'       Correct for one-way light time (also
%                             called "planetary aberration") using a
%                             Newtonian formulation. This correction
%                             yields the location of the terminator
%                             points at the approximate time they
%                             emitted photons arriving at the
%                             observer at 'et' (the difference between
%                             light time to the target center and
%                             light time to the terminator points
%                             is ignored).
%
%                             The light time correction uses an
%                             iterative solution of the light time
%                             equation. The solution invoked by the
%                             'LT' option uses one iteration.
%
%                             The target position as seen by the
%                             observer, the position of the light
%                             source as seen from the target at
%                             et-lt, and the rotation of the target
%                             body, are corrected for light time.
%
%                  'LT+S'     Correct for one-way light time and
%                             stellar aberration using a Newtonian
%                             formulation. This option modifies the
%                             positions obtained with the 'LT' option
%                             to account for the observer's velocity
%                             relative to the solar system
%                             barycenter. This correction also
%                             applies to the position of the light
%                             source relative to the target. The
%                             result is the apparent terminator as
%                              seen by the observer.
%
%                  'CN'       Converged Newtonian light time
%                             correction. In solving the light time
%                             equation, the 'CN' correction iterates
%                             until the solution converges. The
%                             position and rotation of the target
%                             body and the position of the light
%                             source relative to the target are
%                             corrected for light time.
%
%                  'CN+S'     Converged Newtonian light time
%                             and stellar aberration corrections.
%
%               [1,c5] = size(abcorr); char = class(abcorr)
%
%      obsrvr   string name of the observing body. This is typically
%               a spacecraft, the Earth, or a surface point on the
%               Earth. 'obsrvr' is case-insensitive, and leading and
%               trailing blanks in 'obsrvr' are not significant.
%               Optionally, you may supply a string containing the
%               integer ID code for the object. For example both
%               'EARTH' and '399' are legitimate strings that indicate
%               the Earth is the observer.
%
%               [1,c5] = size(obsrvr); char = class(obsrvr)
%
%      npts     number of terminator points to compute.
%
%               [1,1] = size(npts); int32 = class(npts)
%
%   the call:
%
%     [ trgepc, obspos, termpts] = cspice_edterm( trmtyp, source, ...
%                                                 target, et,     ...
%                                                 fixfrm, abcorr, ...
%                                                 obsrvr, npts)
%
%   returns:
%
%      trgepc   the "target epoch" of the calculation. 'trgepc' is
%               defined as follows: letting 'lt' be the one-way light
%               time between the target center and observer, 'trgepc' is
%               either the epoch et-lt or 'et' depending on whether the
%               requested aberration correction is, respectively, for
%               received radiation or omitted. 'lt' is computed using the
%               method indicated by 'abcorr'.
%
%               'trgepc' is expressed as seconds past J2000 TDB.
%
%               [1,1] = size(trgepc); double = class(trgepc)
%
%      obspos   position vector from the center of the target body
%               at epoch 'trgepc' to the observer at epoch 'et'. 'obspos' is
%               expressed in the target body-fixed reference frame
%               'fixref', which is evaluated at 'trgepc'.
%
%               'obspos' is returned to simplify various related
%               computations that would otherwise be cumbersome. For
%               example, the vector 'xvec' from the observer to the
%               ith terminator point can be calculated via the call
%
%                  xvec = trmpts(*,i) - obspos
%
%               To transform the vector 'obspos' from a reference frame
%               'fixref' at time 'trgepc' to a time-dependent reference
%               frame 'ref' at time 'et', the routine pxfrm2_c should be
%               called. Let 'xform' be the 3x3 matrix representing the
%               rotation from the reference frame 'fixref' at time
%               'trgepc' to the reference frame 'ref' at time 'et'. Then
%               'obspos' can be transformed to the result 'refvec' as
%               follows:
%
%                   xform  = cspice_pxfrm2( fixref, ref, trgepc, et )
%                   refvec = xform*obspos
%
%               [3,1] = size(obspos); double = class(obspos)
%
%      trmpts   array of points on the umbral or penumbral terminator
%               of the ellipsoid, as specified by the input argument
%               'trmtyp'. The ith point is contained in the array
%
%                   pos_i = trmpts(*,i)
%
%               Each terminator point is the point of tangency of a
%               plane that is also tangent to the light source. These
%               associated points of tangency on the light source have
%               uniform distribution in longitude when expressed in a
%               cylindrical coordinate system whose Z-axis is the target
%               center to source center vector. The magnitude of the
%               separation in longitude between the tangency points on
%               the light source is
%
%                  2*pi / npts
%
%               If the target is spherical, the terminator points
%               also are uniformly distributed in longitude in the
%               cylindrical system described above. If the target is
%               non-spherical, the longitude distribution of the
%               points generally is not uniform.
%
%               The terminator points are expressed in the body-fixed
%               reference frame designated by 'fixref'. Units are km.
%
%               [3,npts] = size(trmpts); double = class(trmpts)
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
%      Compute sets of umbral and penumbral terminator points on the
%      Moon. Perform a consistency check using the solar incidence
%      angle at each point. We expect to see a solar incidence angle of
%      approximately 90 degrees. Since the solar incidence angle is
%      measured between the local outward normal and the direction to
%      the center of the Sun, the solar incidence angle at an umbral
%      terminator point should exceed 90 degrees by approximately the
%      angular radius of the Sun, while the angle at a penumbral
%      terminator points should be less than 90 degrees by that amount.
%
%      META    = 'standard.tm';
%      NPTS    =  3;
%      first   = true;
%      trmtyps = { 'UMBRAL', 'PENUMBRAL' };
%      s       = [ -1, 1];
%      R2D     = cspice_dpr();
%
%      %
%      % Load meta-kernel.
%      %
%      cspice_furnsh( META )
%
%      %
%      % Set observation time.
%      %
%      utc    = '2007 FEB 3 00:00:00.000';
%
%      et = cspice_str2et( utc );
%
%      %
%      % Set participating objects, frame, and aberration
%      % corrections.
%      %
%      obsrvr = 'EARTH';
%      target = 'MOON';
%      source = 'SUN';
%      fixref = 'IAU_MOON';
%      abcorr = 'LT+S';
%
%      %
%      % Look up the radii of the sun.
%      %
%      srcrad = cspice_bodvrd( source, 'RADII', 3 );
%
%      %
%      % Compute terminator points.
%      %
%      for trmidx=1:2
%
%         [ trgepc, obspos, trmpts] = cspice_edterm(      ...
%                        trmtyps(trmidx), source, target, ...
%                        et,              fixref, abcorr, ...
%                        obsrvr,          NPTS );
%
%         %
%         % Validate terminator points.
%         %
%         % Look up the target-sun vector at the light-time
%         % corrected target epoch.
%         %
%         if ( first )
%            [srcpos, ltime] = cspice_spkpos( source, trgepc, ...
%                                             fixref, abcorr, ...
%                                             target );
%
%            first = false;
%         end
%
%         fprintf(' Terminator type: %s\n', char(trmtyps(trmidx)) )
%
%         for i = 1:NPTS
%
%            %
%            % Convert the ith terminator point to latitudinal
%            % coordinates. Display the point.
%            %
%            [radius, lon, lat] = cspice_reclat( trmpts(:,i) );
%
%            fprintf('Terminator point :%d\n', i )
%            fprintf('  Radius                     (km):  %f\n', radius)
%            fprintf('  Planetocentric longitude   (deg): %f\n', lon *R2D)
%            fprintf('  Planetocentric latitude    (deg): %f\n', lat *R2D)
%
%            %
%            % Find the illumination angles at the
%            % ith terminator point.
%            %
%            [trgepc, srfvec, phase, solar, emissn] = ...
%                                     cspice_ilumin( 'Ellipsoid', ...
%                                            target, et,          ...
%                                            fixref, abcorr,      ...
%                                            obsrvr, trmpts(:,i) );
%
%            fprintf('  Solar incidence angle      (deg): %f\n', solar *R2D)
%
%
%            %
%            % Find the angular radius of the Sun as seen from
%            % the terminator point.
%            %
%            angrad = asin( srcrad(1)/cspice_vdist( srcpos, trmpts(:,i)) );
%
%
%            %
%            % Display the solar incidence angle after
%            % adjusting the angular radius of the Sun
%            % as seen from the terminator point.The
%            % result should be approximately 90 degrees.
%            %
%            fprintf('  Solar incidence angle adjusted for\n' )
%            fprintf('  sun''s angular radius (deg): %18.9f\n\n', ...
%                         ( solar + ( s(trmidx)*angrad ) ) *R2D)
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
%       Terminator type: UMBRAL
%      Terminator point :1
%        Radius                     (km):  1737.400000
%        Planetocentric longitude   (deg): -95.084553
%        Planetocentric latitude    (deg): 0.004053
%        Solar incidence angle      (deg): 90.269766
%        Solar incidence angle adjusted for
%        sun's angular radius (deg):       90.000000125
%
%      Terminator point :2
%        Radius                     (km):  1737.400000
%        Planetocentric longitude   (deg): 84.228092
%        Planetocentric latitude    (deg): 59.995756
%        Solar incidence angle      (deg): 90.269766
%        Solar incidence angle adjusted for
%        sun's angular radius (deg):       90.000000019
%
%      Terminator point :3
%        Radius                     (km):  1737.400000
%        Planetocentric longitude   (deg): 87.216418
%        Planetocentric latitude    (deg): -59.979551
%        Solar incidence angle      (deg): 90.269766
%        Solar incidence angle adjusted for
%        sun's angular radius (deg):       90.000000043
%
%       Terminator type: PENUMBRAL
%      Terminator point :1
%        Radius                     (km):  1737.400000
%        Planetocentric longitude   (deg): 84.914101
%        Planetocentric latitude    (deg): -0.004073
%        Solar incidence angle      (deg): 89.730234
%        Solar incidence angle adjusted for
%        sun's angular radius (deg):       90.000000122
%
%      Terminator point :2
%        Radius                     (km):  1737.400000
%        Planetocentric longitude   (deg): -95.769216
%        Planetocentric latitude    (deg): -59.995785
%        Solar incidence angle      (deg): 89.730234
%        Solar incidence angle adjusted for
%        sun's angular radius (deg):       90.000000021
%
%      Terminator point :3
%        Radius                     (km):  1737.400000
%        Planetocentric longitude   (deg): -92.780892
%        Planetocentric latitude    (deg): 59.979499
%        Solar incidence angle      (deg): 89.730234
%        Solar incidence angle adjusted for
%        sun's angular radius (deg):       90.000000044
%
%-Particulars
%
%   This routine models the boundaries of shadow regions on an
%   ellipsoidal target body "illuminated" by a spherical light
%   source. Light rays are assumed to travel along straight lines;
%   refraction is not modeled.
%
%   Points on the target body's surface are classified according to
%   their illumination as follows:
%
%      -  A target surface point X for which no vector from X to any
%         point in the light source intersects the target, except at
%         X, is considered to be "completely illuminated."
%
%      -  A target surface point X for which each vector from X to a
%         point in the light source intersects the target at points
%         other than X is considered to be "in total shadow."
%
%      -  All other target points are considered to be in partial
%         shadow.
%
%   In this routine, we use the term "umbral terminator" to denote
%   the curve usually called the "terminator": this curve is the
%   boundary of the portion of the target body's surface that lies in
%   total shadow. We use the term "penumbral terminator" to denote
%   the boundary of the completely illuminated portion of the
%   surface.
%
%   In general, the terminator on an ellipsoid is a more complicated
%   curve than the limb (which is always an ellipse). Aside from
%   various special cases, the terminator does not lie in a plane.
%
%   However, the condition for a point X on the ellipsoid to lie on
%   the terminator is simple: a plane tangent to the ellipsoid at X
%   must also be tangent to the light source. If this tangent plane
%   does not intersect the vector from the center of the ellipsoid to
%   the center of the light source, then X lies on the umbral
%   terminator; otherwise X lies on the penumbral terminator.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine edterm_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 18-JUN-2012, EDW (JPL)
%
%-Index_Entries
%
%   find terminator on ellipsoid
%   find umbral terminator on ellipsoid
%   find penumbral terminator on ellipsoid
%
%-&

function [ trgepc, obspos, termpts] = cspice_edterm( trmtyp, source, ...
                                                     target, et,     ...
                                                     fixfrm, abcorr, ...
                                                     obsrvr, npts)

   switch nargin
      case 8

         trmtyp = zzmice_str(trmtyp);
         source = zzmice_str(source);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixfrm = zzmice_str(fixfrm);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         npts   = zzmice_int(npts);

      otherwise

         error ( ['Usage: [ trgepc, obspos, termpts] = cspice_edterm( ' ...
                                       '`trmtyp`, `source`, `target`, ' ...
                                       'et,       `fixfrm`, `abcorr`, ' ...
                                       '`obsrvr`, npts)' ]  )

   end

   %
   % Call the MEX library.
   %
   try
      [ trgepc, obspos, termpts] = mice( 'edterm_c', trmtyp, source, ...
                                                     target, et,     ...
                                                     fixfrm, abcorr, ...
                                                     obsrvr, npts);
   catch
      rethrow(lasterror)
   end


