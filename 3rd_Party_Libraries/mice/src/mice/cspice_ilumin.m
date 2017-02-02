%-Abstract
%
%   CSPICE_ILUMIN computes the illumination angles (phase, solar incidence,
%   and emission) at a specified surface point of a target body.
%
%   This routine supersedes cspice_illum, which lacks an input
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
%      method   parameters defining the computation method to be used.
%               Parameters include, but are not limited to, the shape model
%               used to represent the surface of the target body.
%
%               [1,m] = size(method); char = class(method)
%
%               The only choice currently supported:
%
%                  "Ellipsoid"        The intercept computation uses
%                                     a triaxial ellipsoid to model
%                                     the surface of the target body.
%                                     The ellipsoid's radii must be
%                                     available in the kernel pool.
%
%               Neither case nor white space are significant in
%               `method'. For example, the string ' eLLipsoid ' is
%               valid.
%
%      target   name of the target body. `target' is  case-insensitive, and
%               leading and trailing blanks in `target' are not significant.
%               Optionally, you may supply a string containing the integer ID
%               code for the object. For example both "MOON" and "301" are
%               legitimate strings that indicate the moon is the target body.
%
%               [1,m] = size(target); char = class(target)
%
%      fixref   name of body-fixed, body-centered reference frame associated
%               with the target body. The input surface point `spoint' and the
%               output vector 'srfvec' are expressed relative to this reference
%               frame.
%
%               [1,m] = size(fixref); char = class(fixref)
%
%               'fixref' is case-insensitive, and leading and trailing
%               blanks are not significant.
%
%      et       epoch specified in ephemeris is seconds past J2000, at which
%               the apparent illumination angles at the specified surface point
%               on the target body, as seen from the observing body, are to be
%               computed.
%
%               [1,1] = size(et); double = class(et)
%
%
%      abcorr   aberration correction to be used in computing the location of
%               the surface point, the orientation of the target body, and the
%               location of the Sun.
%
%               [1,m] = size(abcorr); char = class(abcorr)
%
%               For remote sensing applications, where the apparent
%               illumination angles seen by the observer are desired,
%               normally either of the corrections
%
%                  "LT+S"
%                  "CN+S"
%
%               should be used. These and the other supported options
%               are described below. `abcorr' may be any of the
%               following:
%
%                  "NONE"     No aberration correction.
%
%                  "LT"       Correct the position of the input
%                             surface point SPOINT and orientation of
%                             target body for light time, and correct
%                             the position of the Sun as seen
%                             from the target for light time.
%
%                  "LT+S"     Correct the position of the surface
%                             point SPOINT for light time and stellar
%                             aberration, correct the orientation of
%                             the target body for light time, and
%                             correct the position of the Sun as seen
%                             from the target for light time and
%                             stellar aberration.
%
%                  "CN"       Converged Newtonian light time
%                             correction. In solving the light time
%                             equation, the "CN" correction iterates
%                             until the solution converges. Both the
%                             position of the surface point SPOINT c
%                             and rotation of the target body, as
%                             well as the position of the Sun as seen
%                             from the target, are corrected for
%                             light time.
%
%                  "CN+S"     Converged Newtonian light time and
%                             stellar aberration corrections. This
%                             option produces a solution that is at
%                             least as accurate at that obtainable
%                             with the "LT+S" option. Whether the
%                             "CN+S" solution is substantially more
%                             accurate depends on the geometry of the
%                             participating objects and on the
%                             accuracy of the input data. In all
%                             cases this routine will execute more
%                             slowly when a converged solution is
%                             computed.
%
%               Neither case nor white space are significant in
%               `abcorr'. For example, the string
%
%                 "Lt + s"
%
%               is valid.
%
%      obsrvr   name of the observing body. This is typically a  spacecraft,
%               the earth, or a surface point on the earth. `obsrvr' is
%               case-insensitive, and leading and trailing blanks in `obsrvr'
%               are not significant. Optionally, you may supply a string
%               containing the integer ID code for the object. For example both
%               "EARTH" and "399" are legitimate strings that indicate the
%               earth is the observer.
%
%               [1,m] = size(obsrvr); char = class(obsrvr)
%
%               `obsrvr' may be not be identical to `target'.
%
%
%      spoint   a surface location on the target body, expressed in Cartesian
%               coordinates, relative to the body-fixed target frame designated
%               by `fixref'.
%
%               [3,1] = size(spoint); double = class(spoint)
%
%               `spoint' need not be visible from the observer's
%               location at the epoch ET.
%
%               The components of `spoint' have units of km.
%
%   the call:
%
%      [trgepc, srfvec, phase, solar, emissn] = cspice_ilumin( method,  ...
%                                              target, et,     fixref,  ...
%                                              abcorr, obsrvr, spoint)
%
%   returns:
%
%      trgepc   "surface point epoch." 'trgepc' is defined as follows: letting
%               'lt' be the one-way light time between the observer and the
%               input surface point 'spoint', 'trgepc' is either the epoch
%               et-lt or 'et' depending on whether the requested aberration
%               correction is, respectively, for received radiation or omitted.
%               'lt' is computed using the method indicated by 'abcorr'.
%
%               [1,1] = size(trgepc); double = class(trgepc)
%
%               'trgepc' is expressed as seconds past J2000 TDB.
%
%
%      srfvec   the position vector from the observer at 'et' to 'spoint'.
%               'srfvec' is expressed in the target body-fixed  reference frame
%               designated by 'fixref', evaluated at  'trgepc'.
%
%               [3,1] = size(srfvec); double = class(srfvec)
%
%               The components of 'srfvec' are given in units of km.
%
%               One can use the function norm to obtain the
%               distance between the observer and 'spoint':
%
%                     dist = norm( srfvec )
%
%               The observer's position 'obspos', relative to the
%               target body's center, where the center's position is
%               corrected for aberration effects as indicated by
%               'abcorr', can be computed with:
%
%                     obspos = spoint - srfvec
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
%                     xform  = cspice_pxfrm2 ( fixref, ref, trgepc, et )
%                     refvec = xform * srfvec
%
%      phase    phase angle at 'spoint', as seen from 'obsrvr' at time 'et'.
%               This is the angle between the spoint-obsrvr vector and the
%               spoint-sun vector. Units are radians. The range of 'phase' is
%               [0, pi].
%
%               [1,1] = size(phase); double = class(phase)
%
%      solar    solar incidence angle at 'spoint', as seen from 'obsrvr' at
%               time 'et'. This is the angle between the surface normal vector
%               at 'spoint' and the spoint-sun vector. Units are radians. The
%               range of 'solar' is [0, pi].
%
%               [1,1] = size(solar); double = class(solar)
%
%      emissn   emission angle at 'spoint', as seen from 'obsrvr' at time 'et'.
%               This is the angle between the surface normal vector at 'spoint'
%               and the spoint-observer vector. Units are radians. The range of
%               'emissn' is [0, pi].
%
%               [1,1] = size(emissn); double = class(emissn)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Find the phase, solar incidence, and emission angles at the
%      sub-solar and sub-spacecraft points on Mars as seen from the Mars
%      Global Surveyor spacecraft at a user-specified UTC time. Use
%      light time and stellar aberration corrections.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%         KPL/MK
%
%         File name:  standard.tm
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
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0009.tls'
%                                'de421.bsp'
%                                'pck00009.tpc'
%                      )
%
%         \begintext
%
%      %
%      % Load kernel files.
%      %
%      cspice_furnsh( 'standard.tm' )
%      cspice_furnsh( 'mgs_ext13_ipng_mgs95j.bsp' )
%
%      %
%      % Convert the UTC request time to ET (seconds past J2000 TDB).
%      %
%      utc = '2004 JAN 1 12:00:00';
%
%      et = cspice_str2et( utc );
%
%      %
%      % Assign observer and target names. The acronym MGS
%      % indicates Mars Global Surveyor. See NAIF_IDS for a
%      % list of names recognized by SPICE. Also set the
%      % aberration correction flag.
%      %
%      target = 'Mars';
%      obsrvr = 'MGS';
%      abcorr = 'CN+S';
%
%      %
%      % Find the sub-solar point on the Earth as seen from
%      % the MGS spacecraft at et. Use the 'near point'
%      % style of sub-point definition.
%      %
%      [ssolpt, trgepc, srfvec] = ...
%                     cspice_subslr( 'near point: ellipsoid', ...
%                                    target, et, 'iau_mars',  ...
%                                    abcorr, obsrvr );
%
%      %
%      % Now find the sub-spacecraft point.
%      %
%      [sscpt, trgepc, srfvec] = ...
%                     cspice_subpnt( 'near point: ellipsoid', ...
%                                     target, et, 'iau_mars', ...
%                                     abcorr, obsrvr );
%
%      %
%      % Find the phase, solar incidence, and emission
%      % angles at the sub-solar point on the Earth as seen
%      % from MGS at time et.
%      %
%      [ trgepc, srfvec, sslphs, sslsol, sslemi ] = ...
%                     cspice_ilumin( 'Ellipsoid',   ...
%                                     target,  et, 'IAU_MARS', ...
%                                     abcorr,  obsrvr,  ssolpt );
%
%      %
%      % Do the same for the sub-spacecraft point.
%      %
%      [ trgepc, srfvec, sscphs, sscsol, sscemi] = ...
%                      cspice_ilumin( 'Ellipsoid', ...
%                                      target,  et, 'IAU_MARS', ...
%                                      abcorr, obsrvr, sscpt );
%
%      %
%      % Convert the angles to degrees and write them out.
%      %
%      sslphs = sslphs * cspice_dpr;
%      sslsol = sslsol * cspice_dpr;
%      sslemi = sslemi * cspice_dpr;
%      sscphs = sscphs * cspice_dpr;
%      sscsol = sscsol * cspice_dpr;
%      sscemi = sscemi * cspice_dpr;
%
%      fprintf( [ '\n'                                            ...
%                  'UTC epoch is %s\n'                            ...
%                  '\n'                                           ...
%                  'Illumination angles at the sub-solar point:\n'...
%                  '\n'                                           ...
%                  'Phase angle             (deg):  %f\n'         ...
%                  'Solar incidence angle   (deg):  %f\n'         ...
%                  'Emission angle          (deg):  %f\n'         ...
%                  '\n'                                           ...
%                  'The solar incidence angle should be 0.\n'     ...
%                  'The emission and phase angles should be '     ...
%                  'equal.\n'                                     ...
%                  '\n'                                           ...
%                  '\n'                                           ...
%                  'Illumination angles at the sub-s/c point:\n'  ...
%                  '\n'                                           ...
%                  'Phase angle             (deg):  %f\n'         ...
%                  'Solar incidence angle   (deg):  %f\n'         ...
%                  'Emission angle          (deg):  %f\n'         ...
%                  '\n'                                           ...
%                  'The emission angle should be 0.\n'            ...
%                  'The solar incidence and phase angles '        ...
%                  'should be equal.\n'                           ...
%                  '\n' ], ...
%                  utc,    ...
%                  sslphs, ...
%                  sslsol, ...
%                  sslemi, ...
%                  sscphs, ...
%                  sscsol, ...
%                  sscemi );
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      UTC epoch is 2004 JAN 1 12:00:00
%
%      Illumination angles at the sub-solar point:
%
%      Phase angle             (deg):  115.542001
%      Solar incidence angle   (deg):  0.000000
%      Emission angle          (deg):  115.542001
%
%      The solar incidence angle should be 0.
%      The emission and phase angles should be equal.
%
%
%      Illumination angles at the sub-s/c point:
%
%      Phase angle             (deg):  62.084003
%      Solar incidence angle   (deg):  62.084003
%      Emission angle          (deg):  0.000000
%
%      The emission angle should be 0.
%      The solar incidence and phase angles should be equal.
%
%-Particulars
%
%   A sister version of this routine exists named mice_subpnt that returns
%   the output arguments as fields in a single structure.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ilumin_c.
%
%   MICE.REQ
%   FRAMES.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.4, 12-MAR-2012, SCK (JPL), EDW (JPL)
%
%       References to the new 'cspice_pxfrm2' routine were added to the
%       'I/O returns' section. A problem description was added to the
%       'Examples' section.
%
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.3, 29-APR-2011, EDW (JPL)
%
%       Corrected minor typo in header and another in the usage string.
%
%   -Mice Version 1.0.2, 12-MAY-2009, EDW (JPL)
%
%       Edited I/O section; added 'fixref' description.
%
%   -Mice Version 1.0.1, 30-DEC-2008, EDW (JPL)
%
%       Added typography markers to usage string descriptor.
%
%       Minor edit to Example comments.
%
%       Corrected misspellings.
%
%   -Mice Version 1.0.0, 14-FEB-2008, EDW (JPL)
%
%-Index_Entries
%
%   illumination angles
%   lighting angles
%   phase angle
%   emission angle
%   solar incidence angle
%
%-&

function [trgepc, srfvec, phase, solar, emissn] = cspice_ilumin( method, ...
                                                     target, et, fixref, ...
                                                     abcorr, obsrvr, spoint)

   switch nargin
      case 7

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         spoint = zzmice_dp(spoint);

      otherwise

         error ( [ 'Usage: [trgepc, srfvec(3), phase, solar, emissn] = '     ...
                               'cspice_ilumin( `method`, `target`, et, ' ...
                               '`fixref`, `abcorr`, `obsrvr`, spoint(3))' ] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [ilumin] = mice('ilumin_s', method, target, et, ...
                                  fixref, abcorr, obsrvr, spoint);
      trgepc   = reshape( [ilumin(:).trgepc], 1, [] );
      srfvec   = reshape( [ilumin(:).srfvec], 3, [] );
      phase    = reshape( [ilumin(:).phase ], 1, [] );
      solar    = reshape( [ilumin(:).solar ], 1, [] );
      emissn   = reshape( [ilumin(:).emissn], 1, [] );
   catch
      rethrow(lasterror)
   end



