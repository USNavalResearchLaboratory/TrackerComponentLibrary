%-Abstract
%
%   CSPICE_SINCPT computes the surface intercept of the ray on a target
%   body at a specified epoch, optionally corrected for light time and stellar
%   aberration, given an observer and a direction vector defining a ray,
%
%   This routine supersedes cspice_srfxpt, which does not have an input
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
%               the computation method to be used. Parameters
%               include, but are not limited to, the shape model
%               used to represent the surface of the target body.
%
%               The only choice currently supported is
%
%                  "Ellipsoid"        The intercept computation uses
%                                     a triaxial ellipsoid to model
%                                     the surface of the target body.
%                                     The ellipsoid's radii must be
%                                     available in the kernel pool.
%
%               Neither case nor white space are significant in
%               'method'. For example, the string ' eLLipsoid ' is
%               valid.
%
%      target   the scalar string name of the target body. 'target' is
%               case-insensitive, and leading and trailing blanks in
%               'target' are not significant. Optionally, you may
%               supply a string containing the integer ID code
%               for the object. For example both "MOON" and "301"
%               are legitimate strings that indicate the moon is the
%               target body.
%
%
%      et       the double precision scalar epoch of participation of
%               the observer, expressed as ephemeris seconds past J2000
%               TDB: 'et' is the epoch at which the observer's state
%               is computed.
%
%               When aberration corrections are not used, 'et' is also
%               the epoch at which the position and orientation of the
%               target body are computed.
%
%               When aberration corrections are used, 'et' is the epoch
%               at which the observer's state relative to the solar
%               system barycenter is computed; in this case the
%               position and orientation of the target body are
%               computed at et-lt or et+lt, where 'lt' is the one-way
%               light time between the intercept point and the
%               observer, and the sign applied to 'lt' depends on the
%               selected correction. See the description of 'abcorr'
%               below for details.
%
%      fixref   the scalar string name of the body-fixed, body-centered
%               reference frame associated with the target body. The
%               output intercept point 'spoint' and observer to
%               intercept vector 'srfvec' expressed relative to
%               this reference frame.
%
%      abcorr   the scalar string aberration correction to be applied
%               when computing the observer-target state and the
%               orientation of the target body. 'abcorr' may be any of
%               the following.
%
%                  "NONE"     Apply no correction. Return the
%                             geometric surface intercept point on the
%                             target body.
%
%               Let 'lt' represent the one-way light time between the
%               observer and the surface intercept point (note: NOT
%               between the observer and the target body's center).
%               The following values of 'abcorr' apply to the
%               "reception" case in which photons depart from the
%               intercept point's location at the light-time
%               corrected epoch et-lt and *arrive* at the observer's
%               location at 'et':
%
%                  "LT"       Correct for one-way light time (also
%                             called "planetary aberration") using a
%                             Newtonian formulation. This correction
%                             yields the location of the surface
%                             intercept point at the moment it
%                             emitted photons arriving at the
%                             observer at 'et'.
%
%                             The light time correction uses an
%                             iterative solution of the light time
%                             equation. The solution invoked by the
%                             "LT" option uses one iteration.
%
%                             Both the target state as seen by the
%                             observer, and rotation of the target
%                             body, are corrected for light time.
%
%                  "LT+S"     Correct for one-way light time and
%                             stellar aberration using a Newtonian
%                             formulation. This option modifies the
%                             state obtained with the "LT" option to
%                             account for the observer's velocity
%                             relative to the solar system
%                             barycenter. The result is the apparent
%                             surface intercept point as seen by the
%                             observer.
%
%                  "CN"       Converged Newtonian light time
%                             correction. In solving the light time
%                             equation, the "CN" correction iterates
%                             until the solution converges. Both the
%                             state and rotation of the target body
%                             are corrected for light time.
%
%                  "CN+S"     Converged Newtonian light time
%                             and stellar aberration corrections.
%
%               The following values of 'abcorr' apply to the
%               "transmission" case in which photons *depart* from
%               the observer's location at 'et' and arrive at the
%               intercept point at the light-time corrected epoch
%               ET+LT:
%
%                  "XLT"      "Transmission" case: correct for
%                             one-way light time using a Newtonian
%                             formulation. This correction yields the
%                             intercept location at the moment it
%                             receives photons emitted from the
%                             observer's location at 'et'.
%
%                             The light time correction uses an
%                             iterative solution of the light time
%                             equation. The solution invoked by the
%                             "LT" option uses one iteration.
%
%                             Both the target state as seen by the
%                             observer, and rotation of the target
%                             body, are corrected for light time.
%
%                  "XLT+S"    "Transmission" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation  This option modifies the
%                             intercept obtained with the "XLT"
%                             option to account for the observer's
%                             velocity relative to the solar system
%                             barycenter.
%
%                  "XCN"      Converged Newtonian light time
%                             correction. This is the same as XLT
%                             correction but with further iterations
%                             to a converged Newtonian light time
%                             solution.
%
%                  "XCN+S"    "Transmission" case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%               Case and embedded blanks are not significant in 'abcorr'.
%
%      obsrvr   the scalar string name of the observing body. This is
%               typically a spacecraft, the earth, or a surface point on
%               the earth. 'obsrvr' is case-insensitive, and leading and
%               trailing blanks in 'obsrvr' are not significant.
%               Optionally, you may supply a string containing the
%               integer ID code for the object. For example both
%               "MOON" and "301" are legitimate strings that indicate
%               the moon is the observer.
%
%      dref     the scalar string name of the reference frame relative to
%               which the input direction vector is expressed. This may be
%               any frame supported by the SPICE system, including
%               built-in frames (documented in the Frames Required
%               Reading) and frames defined by a loaded frame kernel
%               (FK).
%
%               When 'dref' designates a non-inertial frame, the
%               orientation of the frame is evaluated at an epoch
%               dependent on the frame's center and, if the center is
%               not the observer, on the selected aberration
%               correction. See the description of the direction
%               vector 'dvec' for details.
%
%      dvec     the double precision 3x1 array defining the pointing
%               vector emanating from the observer. The intercept
%               with the target body's surface of the ray
%               defined by the observer and 'dvec' is sought.
%
%               'dvec' is specified relative to the reference frame
%               designated by 'dref'.
%
%               Non-inertial reference frames are treated as follows:
%               if the center of the frame is at the observer's
%               location, the frame is evaluated at 'et'. If the
%               frame's center is located elsewhere, then letting
%               'ltcent' be the one-way light time between the observer
%               and the central body associated with the frame, the
%               orientation of the frame is evaluated at et-ltcent,
%               et+ltcent, or 'et' depending on whether the requested
%               aberration correction is, respectively, for received
%               radiation, transmitted radiation, or is omitted.
%               'ltcent' is computed using the method indicated by
%               'abcorr'.
%
%   the call:
%
%      [ spoint, trgepc, srfvec, found] = cspice_sincpt( method, target, ...
%                                                        et,     fixref, ...
%                                                        abcorr, obsrvr, ...
%                                                        dref,   dvec)
%
%   returns:
%
%      spoint   double precision 3x1 array defining surface intercept
%               point on the target body of the ray defined by the observer
%               and the direction vector. If the ray intersects the target
%               body in multiple points, the selected intersection point is
%               the one closest to the observer. The output argument
%               'found' (see below) indicates whether an intercept was
%               found.
%
%               'spoint' is expressed in Cartesian coordinates,
%               relative to the target body-fixed frame designated by
%               FIXFRM. The body-fixed target frame is evaluated at
%               the intercept epoch 'trgepc' (see description below).
%
%               When light time correction is used, the duration of
%               light travel between 'spoint' to the observer is
%               considered to be the one way light time. When both
%               light time and stellar aberration corrections are
%               used, 'spoint' is selected such that, when 'spoint' is
%               corrected for light time and the vector from the
%               observer to the light-time corrected location of
%               'spoint' is corrected for stellar aberration, the
%               resulting vector is parallel to the ray defined by
%               the observer's location and 'dvec'.
%
%               The components of 'spoint' are given in units of km.
%
%      trgepc   the scalar double precision "intercept epoch."  This is the
%               epoch at which the ray defined by 'obsrvr' and 'dvec'
%               intercepts the target surface at 'spoint'. 'trgepc' is defined
%               as follows: letting 'lt' be the one-way light time between
%               the observer and the intercept point, 'trgepc' is the
%               epoch et-lt, et+lt, or 'et' depending on whether the
%               requested aberration correction is, respectively, for
%               received radiation, transmitted radiation, or
%               omitted. 'lt' is computed using the method indicated by
%               'abcorr'.
%
%               'trgepc' is expressed as seconds past J2000 TDB.
%
%      srfvec   a double precision 3x1 array defining the vector
%               from the observer's position at 'et' to
%               'spoint'. 'srfvec' is expressed in the target body-fixed
%               reference frame designated by 'fixref', evaluated at
%               'trgepc'.
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
%               'abcorr', can be computed via the call:
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
%      found    a scalar logical indicating whether or not the ray
%               intersects the target. If an intersection exists
%               'found' will return as true If the ray misses
%               the target, 'found' will return as false.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      The following program computes surface intercept points on
%      Mars for the boresight and FOV boundary vectors of the MGS MOC
%      narrow angle camera. The intercepts are computed for a single
%      observation epoch. Light time and stellar aberration corrections
%      are used. For simplicity, camera distortion is ignored.
%
%      %
%      %  Local variables
%      %
%      abcorr  = 'CN+S';
%      camera  = 'MGS_MOC_NA';
%      fixref  = 'IAU_MARS';
%      method  = 'Ellipsoid';
%      obsrvr  = 'MGS';
%      target  = 'Mars';
%      utc     = '2003 OCT 13 06:00:00 UTC';
%      NCORNR  = 4;
%
%      %
%      % Load kernel files.
%      %
%      cspice_furnsh( 'standard.tm' );
%      cspice_furnsh( { '/kernels/MGS/ik/moc20.ti',                 ...
%                       '/kernels/MGS/sclk/MGS_SCLKSCET.00061.tsc'  ,...
%                       '/kernels/MGS/spk/mgs_ext12_ipng_mgs95j.bsp',...
%                       '/kernels/MGS/ck/mgs_sc_ext12.bc' } )
%
%      %
%      % Convert the UTC request time to ET (seconds past
%      % J2000, TDB).
%      %
%      et = cspice_str2et( utc );
%
%      %
%      % Get the MGS MOC Narrow angle camera (MGS_MOC_NA)
%      % ID code. Then look up the field of view (FOV)
%      % parameters.
%      %
%      [ camid, found ] = cspice_bodn2c( camera );
%
%      if ( ~found )
%         txt = sprintf( [ 'SPICE(NOTRANSLATION) ' ...
%                         'Could not find ID code for instrument %s.' ], ...
%                          camera);
%         error( txt )
%      end
%
%      %
%      % cspice_getfov will return the name of the camera-fixed frame
%      % in the string 'dref', the camera boresight vector in
%      % the array 'bsight', and the FOV corner vectors in the
%      % array 'bounds'.
%      %
%      [shape, dref, bsight, bounds] = cspice_getfov( camid, NCORNR);
%
%      fprintf (  ['\n' ...
%                  'Surface Intercept Locations for Camera\n'  ...
%                  'FOV Boundary and Boresight Vectors\n'      ...
%                  '\n'                                        ...
%                  '   Instrument:             %s\n'           ...
%                  '   Epoch:                  %s\n'           ...
%                  '   Aberration correction:  %s\n'           ...
%                  '\n'],                                      ...
%                  camera, utc, abcorr )
%
%      for i=1:NCORNR+1
%
%         if( i <= NCORNR )
%            fprintf( 'Corner vector %d\n\n', i)
%            dvec = bounds(:,i);
%         end
%
%         if ( i == (NCORNR + 1) )
%            fprintf( 'Boresight vector\n\n' )
%            dvec = bsight;
%         end
%
%         %
%         % Compute the surface intercept point using
%         % the specified aberration corrections.
%         %
%         [ spoint, trgepc, srfvec, found ] =                   ...
%                        cspice_sincpt( method, target,         ...
%                                       et,     fixref, abcorr, ...
%                                       obsrvr, dref,   dvec );
%         if( found )
%
%            %
%            % Compute range from observer to apparent intercept.
%            %
%            dist = vnorm( srfvec );
%
%            %
%            % Convert rectangular coordinates to planetocentric
%            % latitude and longitude. Convert radians to degrees.
%            %
%            [ radius, lon, lat ] = cspice_reclat( spoint );
%
%            lon = lon * cspice_dpr;
%            lat = lat * cspice_dpr;
%
%            %
%            % Display the results.
%            %
%            fprintf( '  Vector in %s frame = \n', dref )
%            fprintf( '   %18.10e %18.10e %18.10e\n', dvec );
%
%            fprintf( [ '\n'                                              ...
%                       '  Intercept:\n'                                  ...
%                       '\n'                                              ...
%                       '     Radius                   (km)  = %18.10e\n' ...
%                       '     Planetocentric Latitude  (deg) = %18.10e\n' ...
%                       '     Planetocentric Longitude (deg) = %18.10e\n' ...
%                       '     Range                    (km)  = %18.10e\n' ...
%                       '\n' ],                                           ...
%                        radius,  lat,  lon,  dist                          )
%         else
%            disp( 'Intercept not found.' )
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
%      Surface Intercept Locations for Camera
%      FOV Boundary and Boresight Vectors
%
%         Instrument:             MGS_MOC_NA
%         Epoch:                  2003 OCT 13 06:00:00 UTC
%         Aberration correction:  CN+S
%
%      Corner vector 1
%
%        Vector in MGS_MOC_NA frame =
%           1.8571383810e-06  -3.8015622659e-03   9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849411359e+03
%           Planetocentric Latitude  (deg) =  -4.8477481924e+01
%           Planetocentric Longitude (deg) =  -1.2347407905e+02
%           Range                    (km)  =   3.8898310366e+02
%
%      Corner vector 2
%
%        Vector in MGS_MOC_NA frame =
%           1.8571383810e-06   3.8015622659e-03   9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849396987e+03
%           Planetocentric Latitude  (deg) =  -4.8481636340e+01
%           Planetocentric Longitude (deg) =  -1.2339882297e+02
%           Range                    (km)  =   3.8897512129e+02
%
%      Corner vector 3
%
%        Vector in MGS_MOC_NA frame =
%          -1.8571383810e-06   3.8015622659e-03   9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849396899e+03
%           Planetocentric Latitude  (deg) =  -4.8481661910e+01
%           Planetocentric Longitude (deg) =  -1.2339882618e+02
%           Range                    (km)  =   3.8897466238e+02
%
%      Corner vector 4
%
%        Vector in MGS_MOC_NA frame =
%          -1.8571383810e-06  -3.8015622659e-03   9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849411271e+03
%           Planetocentric Latitude  (deg) =  -4.8477507498e+01
%           Planetocentric Longitude (deg) =  -1.2347408220e+02
%           Range                    (km)  =   3.8898264472e+02
%
%      Boresight vector
%
%        Vector in MGS_MOC_NA frame =
%           0.0000000000e+00   0.0000000000e+00   1.0000000000e+00
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849404102e+03
%           Planetocentric Latitude  (deg) =  -4.8479579822e+01
%           Planetocentric Longitude (deg) =  -1.2343645396e+02
%           Range                    (km)  =   3.8897573572e+02
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine sincpt_c.
%
%   MICE.REQ
%   FRAMES.REQ
%   NAIF_IDS.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.3, 12-MAR-2012, SCK (JPL)
%
%      References to the new 'cspice_pxfrm2' routine were
%      added to the 'I/O returns' section. A problem description was
%      added to the 'Examples' section, and the references to
%      'srfxpt_c' and the second example were removed.
%
%   -Mice Version 1.0.2, 14-JUL-2010, EDW (JPL)
%
%      Corrected minor typo in header.
%
%   -Mice Version 1.0.1, 23-FEB-2009, EDW (JPL)
%
%      Added proper markers for usage string variable types.
%
%   -Mice Version 1.0.0, 11-FEB-2008, EDW (JPL)
%
%-Index_Entries
%
%   find surface intercept point
%   find intersection of ray and target body surface
%   find intercept of ray on target body surface
%
%-&

function [ spoint, trgepc, srfvec, found] = ...
          cspice_sincpt( method, target, et, fixref, abcorr, obsrvr, dref, dvec)

   switch nargin
      case 8

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         dref   = zzmice_str(dref);
         dvec   = zzmice_dp(dvec);

      otherwise

         error( [ 'Usage: [ spoint, trgepc, srfvec, found] =  ' ...
                  'cspice_sincpt( `method`, `target`, et, `fixref`, ' ...
                                 '`abcorr`, `obsrvr`, `dref`, dvec)' ]  )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [sincpt] = mice('sincpt_s', method, target, ...
                                  et, fixref, abcorr, obsrvr, dref, dvec);
      spoint = reshape( [sincpt.spoint], 3, [] );
      trgepc = reshape( [sincpt.trgepc], 1, [] );
      srfvec = reshape( [sincpt.srfvec], 3, [] );
      found  = reshape( [sincpt.found] , 1, [] );
   catch
      rethrow(lasterror)
   end

