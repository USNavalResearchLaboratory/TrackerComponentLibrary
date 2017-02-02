%-Abstract
%
%   CSPICE_SUBPNT computes the rectangular coordinates of the
%   sub-observer point on a target body at a specified epoch,
%   optionally corrected for light time and stellar aberration.
%
%   This routine supersedes cspice_subpt, which does not have an input
%   argument for the target body-fixed frame name.
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
%      method   a scalar string providing parameters defining
%               the computation method to be used.
%
%               The supported values of 'method' are listed below.
%               Please note that the colon is a required delimiter;
%               using a blank will not work.
%
%                     'Near point: ellipsoid'   The sub-observer point
%                                               computation uses a
%                                               triaxial ellipsoid to
%                                               model the surface of the
%                                               target body. The
%                                               sub-observer point is
%                                               defined as the nearest
%                                               point on the target
%                                               relative to the observer.
%
%                     'Intercept: ellipsoid'    The sub-observer point
%                                               computation uses a
%                                               triaxial ellipsoid to
%                                               model the surface of the
%                                               target body. The
%                                               sub-observer point is
%                                               defined as the target
%                                               surface intercept of the
%                                               line containing the
%                                               observer and the
%                                               target's center.
%
%               Neither case nor white space are significant in
%               'method'. For example, the string
%
%                    ' nearpoint:ELLIPSOID '
%
%               is valid.
%
%      target   the scalar string name of the target body. The target
%               body is an ephemeris object (its trajectory is given by
%               SPK data), and is an extended object.
%
%               The string 'target' is case-insensitive, and leading
%               and trailing blanks in 'target' are not significant.
%               Optionally, you may supply a string containing the
%               integer ID code for the object. For example both
%               'MOON' and '301' are legitimate strings that indicate
%               the moon is the target body.
%
%               When the target body's surface is represented by a
%               tri-axial ellipsoid, this routine assumes that a
%               kernel variable representing the ellipsoid's radii is
%               present in the kernel pool. Normally the kernel
%               variable would be defined by loading a PCK file.
%
%
%      et       the double precision scalar epoch, expressed as seconds
%               past J2000 TDB, of the observer: 'et' is
%               the epoch at which the observer's state is computed.
%
%               When aberration corrections are not used, 'et' is also
%               the epoch at which the position and orientation of
%               the target body are computed.
%
%               When aberration corrections are used, 'et' is the epoch
%               at which the observer's state relative to the solar
%               system barycenter is computed; in this case the
%               position and orientation of the target body are
%               computed at et-lt or et+lt, where 'lt' is the one-way
%               light time between the sub-observer point and the
%               observer, and the sign applied to 'lt' depends on the
%               selected correction. See the description of 'abcorr'
%               below for details.
%
%      fixref   the scalar string name of the body-fixed, body-centered
%               reference frame associated with the target body.
%               The output sub-observer point 'spoint' will be
%               expressed relative to this reference frame.
%
%      abcorr   the scalar string aberration correction to apply
%               when computing the observer-target state and the
%               orientation of the target body.
%
%               For remote sensing applications, where the apparent
%               sub-observer point seen by the observer is desired,
%               normally either of the corrections
%
%                     'LT+S'
%                     'CN+S'
%
%               should be used. These and the other supported options
%               are described below. 'abcorr' may be any of the
%               following:
%
%                     'NONE'     Apply no correction. Return the
%                                geometric sub-observer point on the
%                                target body.
%
%               Let 'lt' represent the one-way light time between the
%               observer and the sub-observer point (note: NOT
%               between the observer and the target body's center).
%               The following values of 'abcorr' apply to the
%               "reception" case in which photons depart from the
%               sub-observer point's location at the light-time
%               corrected epoch et-lt and *arrive* at the observer's
%               location at 'et':
%
%                     'LT'       Correct for one-way light time (also
%                                called "planetary aberration") using a
%                                Newtonian formulation. This correction
%                                yields the location of sub-observer
%                                point at the moment it emitted photons
%                                arriving at the observer at 'et'.
%
%                                The light time correction uses an
%                                iterative solution of the light time
%                                equation. The solution invoked by the
%                                'LT' option uses one iteration.
%
%                                Both the target position as seen by the
%                                observer, and rotation of the target
%                                body, are corrected for light time.
%
%                     'LT+S'     Correct for one-way light time and
%                                stellar aberration using a Newtonian
%                                formulation. This option modifies the
%                                state obtained with the 'LT' option to
%                                account for the observer's velocity
%                                relative to the solar system
%                                barycenter. The result is the apparent
%                                sub-observer point as seen by the
%                                observer.
%
%                     'CN'       Converged Newtonian light time
%                                correction. In solving the light time
%                                equation, the 'CN' correction iterates
%                                until the solution converges. Both the
%                                position and rotation of the target
%                                body are corrected for light time.
%
%                     'CN+S'     Converged Newtonian light time and
%                                stellar aberration corrections. This
%                                option produces a solution that is at
%                                least as accurate at that obtainable
%                                with the 'LT+S' option. Whether the 'CN+S'
%                                solution is substantially more accurate
%                                depends on the geometry of the
%                                participating objects and on the
%                                accuracy of the input data. In all
%                                cases this routine will execute more
%                                slowly when a converged solution is
%                                computed.
%
%               The following values of 'abcorr' apply to the
%               "transmission" case in which photons *depart* from
%               the observer's location at 'et' and arrive at the
%               sub-observer point at the light-time corrected epoch
%               et+lt:
%
%                     'XLT'      "Transmission" case: correct for
%                                 one-way light time using a Newtonian
%                                 formulation. This correction yields the
%                                sub-observer location at the moment it
%                                receives photons emitted from the
%                                observer's location at 'et'.
%
%                                The light time correction uses an
%                                iterative solution of the light time
%                                equation. The solution invoked by the
%                                'LT' option uses one iteration.
%
%                                Both the target position as seen by the
%                                observer, and rotation of the target
%                                body, are corrected for light time.
%
%                     'XLT+S'    "Transmission" case: correct for
%                                one-way light time and stellar
%                                aberration using a Newtonian
%                                formulation  This option modifies the
%                                sub-observer point obtained with the
%                                'XLT' option to account for the
%                                observer's velocity relative to the
%                                solar system barycenter.
%
%                     'XCN'      Converged Newtonian light time
%                                correction. This is the same as XLT
%                                correction but with further iterations
%                                to a converged Newtonian light time
%                                solution.
%
%                     'XCN+S'    "Transmission" case: converged
%                                Newtonian light time and stellar
%                                aberration corrections.
%
%      obsrvr   the scalar string name of the observing body. The
%               observing body is an ephemeris object: it typically
%               is a spacecraft, the earth, or a surface point on the
%               earth. 'obsrvr' is case-insensitive, and leading and
%               'obsrvr' are not significant. Optionally, you may
%               trailing blanks in supply a string containing the integer
%               ID code for the object. For example both 'MOON' and '301'
%               are legitimate strings that indicate the Moon is the
%               observer.
%
%   the call:
%
%      [spoint, trgepc, srfvec] = cspice_subpnt( method, target, ...
%                                                et,     fixref, ...
%                                                abcorr, obsrvr)
%
%   returns:
%
%      spoint   a double precision 3x1 array defining the sub-observer point
%               on the target body.
%
%               The sub-observer point is defined either as the point
%               on the target body that is closest to the observer,
%               or the target surface intercept of the line from the
%               observer to the target's center; the input argument
%               'method' selects the definition to be used.
%
%               'spoint' is expressed in Cartesian coordinates,
%               relative to the body-fixed target frame designated by
%               'fixref'. The body-fixed target frame is evaluated at
%               the sub-observer epoch 'trgepc' (see description below).
%
%               When light time correction is used, the duration of
%               light travel between 'spoint' to the observer is
%               considered to be the one way light time.
%
%               When aberration corrections are used, 'spoint' is
%               computed using target body position and orientation
%               that have been adjusted for the corrections
%               applicable to 'spoint' itself rather than to the target
%               body's center. In particular, if the stellar
%               aberration correction applicable to 'spoint' is
%               represented by a shift vector 's', then the light-time
%               corrected position of the target is shifted by 's'
%               before the sub-observer point is computed.
%
%               The components of 'spoint' have units of km.
%
%      trgepc   the scalar double precision "sub-observer point epoch."
%              'trgepc' is defined as follows: letting 'lt' be the one-way
%               light time between the observer and the sub-observer point,
%               'trgepc' is the epoch et-lt, et+lt, or 'et' depending on
%               whether the requested aberration correction is,
%               respectively, for received radiation, transmitted
%               radiation, or omitted. 'lt' is computed using the
%               method indicated by 'abcorr'.
%
%              'trgepc' is expressed as seconds past J2000 TDB.
%
%      srfvec   a double precision 3x1 array defining the position vector
%               from the observer at 'et' to 'spoint'. 'srfvec'
%               is expressed in the target body-fixed  reference frame
%               designated by 'fixref', evaluated at  'trgepc'.
%
%               The components of 'srfvec' are given in units of km.
%
%               One can use the CSPICE function vnorm_c to obtain the
%               distance between the observer and 'spoint':
%
%                  dist = norm( srfvec )
%
%               The observer's position 'obspos', relative to the
%               target body's center, where the center's position is
%               corrected for aberration effects as indicated by
%               'abcorr', can be computed with:
%
%                  obspos = spoint - srfvec
%
%               To transform the vector 'srfvec' from a reference frame
%               'fixref' at time 'trgepc' to a time-dependent reference
%               frame 'ref' at time 'et', the routine 'cspice_pxfrm2' should be
%               called. Let 'xform' be the 3x3 matrix representing the
%               rotation from the reference frame 'fixref' at time
%               'trgepc' to the reference frame 'ref' at time 'et'. Then
%               'srfvec' can be transformed to the result 'refvec' as
%               follows:
%
%                  xform  = cspice_pxfrm2 ( fixref, ref, trgepc, et )
%                  refvec = xform * srfvec
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Find the sub-Earth point on Mars for a specified time. Perform
%      the computation twice, using both the "intercept" and "near
%      point" options. Display the location of both the Earth and the
%      sub-Earth point using both planetocentric and planetographic
%      coordinates.
%
%      %
%      % Load kernel files via the meta-kernel.
%      %
%      cspice_furnsh( 'standard.tm' );
%
%      %
%      % Convert the UTC request time to ET (seconds past
%      % J2000, TDB).
%      %
%      et = cspice_str2et( '2008 aug 11 00:00:00' );
%
%      %
%      % Look up the target body's radii. We'll use these to
%      % convert Cartesian to planetodetic coordinates. Use
%      % the radii to compute the flattening coefficient of
%      % the reference ellipsoid.
%      %
%      radii  = cspice_bodvrd( 'MARS', 'RADII', 3 );
%
%      %
%      % Let RE and RP be, respectively, the equatorial and
%      % polar radii of the target.
%      %
%      re = radii(1);
%      rp = radii(3);
%      f = ( re-rp)/re;
%
%      %
%      % Compute sub-observer point using light time and stellar
%      % aberration corrections. Use the "target surface intercept"
%      % definition of the sub-observer point on the first loop
%      % iteration, and use the "near point" definition on the
%      % second.
%      %
%
%      method = { 'Intercept:  ellipsoid', 'Near point: ellipsoid' };
%
%      for i=1:2
%
%         [spoint, trgepc, srfvec] = cspice_subpnt( method(i), ...
%                         'MARS', et, 'IAU_MARS', 'LT+S', 'EARTH' );
%
%         %
%         % Compute the observer's distance from SPOINT.
%         %
%         odist = norm(srfvec);
%
%         %
%         % Convert the sub-observer point's rectangular coordinates
%         % to planetographic longitude, latitude and altitude.
%         % Convert radians to degrees.
%         %
%         [ spglon, spglat, spgalt] = cspice_recpgr( 'mars', spoint, re, f );
%
%         spglon = spglon * cspice_dpr;
%         spglat = spglat * cspice_dpr;
%
%         %
%         % Convert sub-observer point's rectangular coordinates to
%         % planetocentric radius, longitude, and latitude. Convert
%         % radians to degrees.
%         %
%         [ spcrad, spclon, spclat ] =cspice_reclat( spoint ) ;
%
%         spclon = spclon * cspice_dpr;
%         spclat = spclat * cspice_dpr;
%
%         %
%         % Compute the observer's position relative to the center of the
%         % target, where the center's location has been adjusted using
%         % the aberration corrections applicable to the sub-point.
%         % Express the observer's location in geodetic coordinates.
%         %
%         obspos = spoint - srfvec;
%
%         [ opglon, opglat, opgalt] = cspice_recpgr( 'mars', obspos, re, f );
%
%         opglon = opglon * cspice_dpr;
%         opglat = opglat * cspice_dpr;
%
%         %
%         % Convert the observer's rectangular coordinates to planetodetic
%         % longitude, latitude and altitude. Convert radians to degrees.
%         %
%         [opcrad, opclon, opclat] = cspice_reclat( obspos ) ;
%
%         opclon = opclon * cspice_dpr;
%         opclat = opclat * cspice_dpr;
%
%         fprintf( 'Computational Method %s\n\n', char(method(i)) )
%
%         fprintf( 'Observer altitude                      (km) = %21.9f\n', ...
%                                                                   opgalt )
%         fprintf( 'Length of SRFVEC                       (km) = %21.9f\n', ...
%                                                                   odist  )
%         fprintf( 'Sub-observer point altitude            (km) = %21.9f\n', ...
%                                                                   spgalt )
%         fprintf( 'Sub-observer planetographic longitude (deg) = %21.9f\n', ...
%                                                                   spglon )
%         fprintf( 'Observer planetographic longitude     (deg) = %21.9f\n', ...
%                                                                   opglon )
%         fprintf( 'Sub-observer planetographic latitude  (deg) = %21.9f\n', ...
%                                                                   spglat )
%         fprintf( 'Observer planetographic latitude      (deg) = %21.9f\n', ...
%                                                                   opglat )
%         fprintf( 'Sub-observer planetocentric longitude (deg) = %21.9f\n', ...
%                                                                   spclon )
%         fprintf( 'Observer planetocentric longitude     (deg) = %21.9f\n', ...
%                                                                   opclon )
%         fprintf( 'Sub-observer planetocentric latitude  (deg) = %21.9f\n', ...
%                                                                   spclat )
%         fprintf( 'Observer planetocentric latitude      (deg) = %21.9f\n', ...
%                                                                   opclat )
%         fprintf( '\n')
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
%      Computational Method Intercept:  ellipsoid
%
%      Observer altitude                      (km) =   349199089.604656994
%      Length of SRFVEC                       (km) =   349199089.641352594
%      Sub-observer point altitude            (km) =          -0.000000000
%      Sub-observer planetographic longitude (deg) =         199.302304818
%      Observer planetographic longitude     (deg) =         199.302304818
%      Sub-observer planetographic latitude  (deg) =          26.262401078
%      Observer planetographic latitude      (deg) =          25.994936593
%      Sub-observer planetocentric longitude (deg) =         160.697695182
%      Observer planetocentric longitude     (deg) =         160.697695182
%      Sub-observer planetocentric latitude  (deg) =          25.994934013
%      Observer planetocentric latitude      (deg) =          25.994934013
%
%      Computational Method Near point: ellipsoid
%
%      Observer altitude                      (km) =   349199089.604648590
%      Length of SRFVEC                       (km) =   349199089.604648590
%      Sub-observer point altitude            (km) =          -0.000000000
%      Sub-observer planetographic longitude (deg) =         199.302304819
%      Observer planetographic longitude     (deg) =         199.302304819
%      Sub-observer planetographic latitude  (deg) =          25.994936593
%      Observer planetographic latitude      (deg) =          25.994936593
%      Sub-observer planetocentric longitude (deg) =         160.697695181
%      Observer planetocentric longitude     (deg) =         160.697695181
%      Sub-observer planetocentric latitude  (deg) =          25.729407071
%      Observer planetocentric latitude      (deg) =          25.994934013
%
%-Particulars
%
%   A sister version of this routine exists named mice_subpnt that returns
%   the output arguments as fields in a single structure.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine subpnt_c.
%
%   MICE.REQ
%   FRAMES.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.2, 25-OCT-2011, SCK (JPL)
%
%      References to the new 'cspice_pxfrm2' routine were
%      added to the 'I/O returns' section.
%
%   -Mice Version 1.0.1, 12-MAY-2009, EDW (JPL)
%
%      Corrected type in I/O call description. The call description
%      lacked the 'fixref' argument.
%
%   -Mice Version 1.0.0, 30-JAN-2008, EDW (JPL)
%
%-Index_Entries
%
%   find sub-observer point on target body
%   find nearest point to observer on target body
%
%-&

function [spoint, trgepc, srfvec] = cspice_subpnt( method, target, et, ...
                                                   fixref, abcorr, obsrvr )

   switch nargin
      case 6

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);

      otherwise

         error ( ['Usage: [spoint, trgepc, srfvec] = '     ...
                  'cspice_subpnt( `method`, `target`,'     ...
                  ' et, `fixref`, `abcorr`, `obsrvr`)']  )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [subpnt] = mice('subpnt_s', method, target, et, fixref, abcorr, obsrvr);
      spoint   = reshape( [subpnt.spoint], 3, [] );
      trgepc   = reshape( [subpnt.trgepc], 1, [] );
      srfvec   = reshape( [subpnt.srfvec], 3, [] );
   catch
      rethrow(lasterror)
   end



