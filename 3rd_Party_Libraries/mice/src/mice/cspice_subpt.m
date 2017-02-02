%-Abstract
%
%   CSPICE_SUBPT determines the coordinates of the sub-observer point
%   on a target body at a particular epoch, optionally corrected
%   for planetary (light time) and stellar aberration. The call also
%   returns the observer's altitude above the target body.
%
%   Deprecated: This routine has been superseded by the routine
%   cspice_subpnt. This routine is supported for purposes of
%   backward compatibility only.
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
%      method   the scalar string name of the computation method.
%               The choices are:
%
%                  'Near point'       The sub-observer point is
%                                     defined as the nearest point on
%                                     the target relative to the
%                                     observer.
%
%                  'Intercept'        The sub-observer point is
%                                     defined as the target surface
%                                     intercept of the line
%                                     containing the observer and the
%                                     target's center.
%
%               In both cases, the intercept computation treats the
%               surface of the target body as a triaxial ellipsoid.
%               The ellipsoid's radii must be available in the kernel
%               pool.
%
%               Neither case nor white space are significant in
%               'method'.  For example, the string ' NEARPOINT' is
%               valid.
%
%      target   the scalar string name of the observed target body. 'target'
%               is case-insensitive, and leading and trailing blanks in
%               'target' are not significant. Optionally, you may supply
%               a string containing the integer ID code for the object.
%               For example both 'MOON' and '301' are legitimate strings
%               that indicate the moon is the target body.
%
%               This routine assumes that the target body is modeled by
%               a tri-axial ellipsoid, and that a PCK file containing
%               its radii has been loaded into the kernel pool via
%               cspice_furnsh.
%
%      et       the double precision scalar or 1xN array of epochs in
%               ephemeris seconds past J2000 at which to compute the
%               sub-observer point on the target body.
%
%      abcorr   the scalar string indicating the aberration corrections to
%               apply when computing the observer-target state.  'abcorr'
%               may be any of the following.
%
%                  'NONE'     Apply no correction. Return the
%                             geometric sub-observer point on the
%                             target body.
%
%                  'LT'       Correct for planetary (light time)
%                             aberration.  Both the state and rotation
%                             of the target body are corrected for
%                             light time.
%
%                  'LT+S'     Correct for planetary (light time) and
%                             stellar aberration. Both the state and
%                             rotation of the target body are
%                             corrected for light time.
%
%                  'CN'       Converged Newtonian light time correction.
%                             This option produces a solution that is at
%                             least as accurate at that obtainable 
%                             with the 'LT' option. Whether the 'CN' 
%                             solution is substantially more accurate 
%                             depends on the geometry of the 
%                             participating objects and on the 
%                             accuracy of the input data. In all 
%                             cases this routine will execute more 
%                             slowly when a converged solution is 
%                             computed. See the section titled "The 
%                             Computation of Light Time" in the SPK
%                             Required Reading document spk.req for 
%                             details.
%
%                  'CN+S'     Converged Newtonian light time
%                             and stellar aberration corrections. 
%                             Both the state and rotation of the
%                             target body are corrected for light
%                             time.
%
%      obsrvr   the scalar string name of the observing body.  This is
%               typically a spacecraft, the earth, or a surface point on the
%               earth. 'obsrvr' is case-insensitive, and leading and trailing
%               blanks in 'obsrvr' are not significant. Optionally, you may
%               supply a string containing the integer ID code for the object.
%               For example both 'EARTH' and '399' are legitimate strings that
%               indicate the earth is the observer.
%
%   the call:
%
%      [spoint, alt] = cspice_subpt( method, target, et, abcorr, obsrvr)
%
%   returns:
%
%      spoint   a double precision 3x1 array or 3xN array of
%               coordinates of the 'obsrvr' subpoint on 'target'
%               at 'et'  expressed relative to the body-fixed frame of the
%               target body.
%
%               The sub-observer point is defined either as the point
%               on the target body that is closest to the observer,
%               or the target surface intercept of the line from the
%               observer to the target's center; the input argument
%               'method' selects the definition to be used.
%
%               The body-fixed frame, which is time-dependent, is
%               evaluated at 'et' if 'abcorr' is 'NONE'; otherwise the
%               frame is evaluated at et-lt, where 'lt' is the one-way
%               light time from target to observer.
%
%               The state of the target body is corrected for
%               aberration as specified by 'abcorr'; the corrected
%               state is used in the geometric computation.  As
%               indicated above, the rotation of the target is
%               retarded by one-way light time if 'abcorr' specifies
%               that light time correction is to be done.
%
%      alt      the double precision scalar or Nx1 array of altitudes of
%               'obsrvr' above 'target'. When 'method' specifies a
%               "near point" computation, 'alt' is truly altitude in the
%               standard geometric sense:  the length of a segment dropped
%               from the observer to the target's surface, such that the
%               segment is perpendicular to the surface at the
%               contact point 'spoint'.
%
%               When 'method' specifies an "intercept" computation, 'alt'
%               is still the length of the segment from the observer
%               to the surface point 'spoint', but this segment in
%               general is not perpendicular to the surface.
%               above 'target'
%
%               'spoint' and 'alt' return with the same vectorization
%                measure (N) as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      %
%      % Find the sub point position of the moon on the earth at
%      % epoch JAN 1, 2006 using the "near point" then the "intercept"
%      % options. Apply light time correction to return apparent position.
%      %
%      %
%      % Load the meta kernel listing the needed SPK, PCK, LSK
%      % kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Calculate the location of the sub point of the moon
%      % on the earth .
%
%      et = cspice_str2et( 'JAN 1, 2006' );
%
%      %
%      % First use option 'Near Point'
%      %
%      [point1,alt1] = cspice_subpt( 'near point', 'earth', et, 'lt+s', 'moon');
%
%      disp( 'Sub-point location  coordinates - near point:' )
%      fprintf( '    %15.8f\n', point1 )
%
%      disp( 'Sub-point observer altitude:' )
%      fprintf( '    %15.8f\n', alt1 )
%
%      disp(' ')
%
%      %
%      % Now use option 'Intercept'
%      %
%      [point2,alt2] = cspice_subpt( 'intercept', 'earth', et, 'lt+s', 'moon');
%
%      disp( 'Sub-point location coordinates - intercept:' )
%      fprintf( '    %15.8f\n', point2 )
%
%      disp( 'Sub-point observer altitude:' )
%      fprintf( '    %15.8f\n', alt2 )
%
%      %
%      % Calculate the Euclidean distance between the two locations
%      % and the angular separation between the position vectors.
%      %
%      dist = norm( point1 - point2);
%      sep  = cspice_vsep(point1, point2 )*cspice_dpr;
%
%      disp(' ')
%
%      fprintf( 'Distance between locations            (km): %8.5f\n', dist);
%      fprintf( 'Angular separation between locations (deg): %8.5f\n', sep );
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%        Sub-point location  coordinates - near point:
%             -5532.84463404
%             -1443.48660124
%             -2816.23526241
%        Sub-point observer altitude:
%            356091.70776573
%
%        Sub-point location coordinates - intercept
%             -5525.64311958
%             -1441.60776851
%             -2831.19585471
%        Sub-point observer altitude:
%            356091.73014965
%
%        Distance between locations            (km): 16.70961
%        Angular separation between locations (deg):  0.15020
%
%   Example(2):
%
%      %
%      % Find the sub body position of the moon on the earth at
%      % at epoch JAN 1, 2006 and for the next 12 months. Use the
%      % 'near point' option to calculate the physically
%      % closest point between the two bodies.
%      %
%      % Load the meta kernel listing the needed SPK, PCK, LSK
%      % kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Convert the calendar string to ephemeris time.
%      %
%      et0 = cspice_str2et( 'JAN 1, 2006' );
%
%      %
%      % Fill an array with epochs, start with the epoch above then
%      % epochs in steps on one month ( thirty days in seconds)
%      %
%      et  = [0:12]*cspice_spd*30. + et0;
%
%      %
%      % Calculate the nearpoint of the moon with respect to earth at
%      % the epochs defined in 'et'.
%      %
%      [point,alt] = cspice_subpt( 'near point', 'earth', et, 'lt+s', 'moon');
%
%      %
%      % Convert the subpoint coordinates to lat/lon expressed in degrees with
%      % the radius.
%      %
%      % Note, 'radius' and 'alt' do not refer to the same quantity.
%      %
%      [radius, longitude, latitude] = cspice_reclat(point);
%      longitude                     = longitude * cspice_dpr;
%      latitude                      = latitude  * cspice_dpr;
%
%      %
%      % Convert the 'et' epochs to calendar format.
%      %
%      utc = cspice_et2utc( et, 'C', 3 );
%
%      for n=1:13
%         txt = sprintf( 'Moon subpoint epoch: %s', utc(n,:) );
%         disp( txt )
%
%         txt = sprintf( '              (deg): longitude %8.4f', longitude(n) );
%         disp( txt )
%
%         txt = sprintf( '              (deg): latitude  %8.4f', latitude(n) );
%         disp( txt )
%         disp( ' ' )
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
%                 ... partial output ...
%
%      Moon subpoint epoch: 2006 JUL 30 00:00:00.001
%                    (deg): longitude -127.7548
%                    (deg): latitude   -0.1948
%
%      Moon subpoint epoch: 2006 AUG 29 00:00:00.001
%                    (deg): longitude -128.2727
%                    (deg): latitude  -15.0349
%
%      Moon subpoint epoch: 2006 SEP 28 00:00:00.002
%                    (deg): longitude -123.9021
%                    (deg): latitude  -25.9738
%
%      Moon subpoint epoch: 2006 OCT 28 00:00:00.001
%                    (deg): longitude -113.7475
%                    (deg): latitude  -27.7753
%
%      Moon subpoint epoch: 2006 NOV 27 00:00:00.001
%                    (deg): longitude -104.0459
%                    (deg): latitude  -17.9194
%
%      Moon subpoint epoch: 2006 DEC 27 00:00:00.000
%                    (deg): longitude -98.2728
%                    (deg): latitude   -0.5411
%
%-Particulars
%
%   A sister version of this routine exists named mice_subpt that returns
%   the output arguments as fields in a single structure.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine subpt_c.
%
%   MICE.REQ
%   FRAMES.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.3, 23-JUN-2014, NJB (JPL)
%
%      Updated description of converged Newtonian light time
%      correction. Replaced double quotes with single quotes
%      in constant strings appearing in comments.
%
%   -Mice Version 1.0.2, 18-MAY-2010, BVS (JPL)
%
%      Index line now states that this routine is deprecated.
%
%   -Mice Version 1.0.1, 11-NOV-2008, EDW (JPL)
%
%      Edits to header; Abstract now states that this routine is
%      deprecated.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   DEPRECATED sub-observer point
%
%-&

function [spoint, alt] = cspice_subpt( method, target, et, abcorr, obsrvr )

   switch nargin
      case 5

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);

      otherwise

         error ( ['Usage: [_spoint(3)_, _alt_] = '     ...
                  'cspice_subpt( `method`, `target`, ' ...
                  '_et_, `abcorr`, `obsrvr`)']  )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [subpt] = mice('subpt_s', method, target, et, abcorr, obsrvr);
      spoint   = reshape( [subpt.pos], 3, [] );
      alt      = reshape( [subpt.alt], 1, [] );
   catch
      rethrow(lasterror)
   end



