%-Abstract
%
%   CSPICE_SUBSOL determines the coordinates of the sub-solar
%   point on a target  body as seen by a specified observer at a
%   specified epoch, optionally corrected for planetary (light time)
%   and stellar aberration.
%
%   Deprecated: This routine has been superseded by the routine
%   cspice_subslr. This routine is supported for purposes of
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
%      method   the scalar sting specifying the computation method
%               to be used.  The choices are:
%
%                  'Near point'    The sub-solar point is defined
%                                  as the nearest point on the
%                                  target to the sun.
%
%                  'Intercept'     The sub-observer point is defined
%                                  as the target surface intercept of
%                                  the line containing the target's
%                                  center and the sun's center.
%
%               In both cases, the intercept computation treats the
%               surface of the target body as a triaxial ellipsoid.
%               The ellipsoid's radii must be available in the kernel
%               pool.
%
%               Neither case nor white space are significant in
%               method.  For example, the string ' NEARPOINT' is
%               valid.
%
%
%      target   the scalar sting describing the name of the target body.
%               'target' is case-insensitive, and leading and trailing
%               blanks in 'target' are not significant. Optionally, you
%               may supply a string containing the integer ID code for
%               the object. For example both 'MOON' and '301' are legitimate
%               strings that indicate the moon is the target body.
%
%               This routine assumes that the target body is modeled by
%               a tri-axial ellipsoid, and that a PCK file containing
%               its radii has been loaded into the kernel pool via
%               cspice_furnsh.
%
%      et       the double precision scalar or 1xN array of ephemeris
%               time expressed as ephemeris seconds past J2000 at which
%               the sub-solar point on the target body is to be
%               computed.
%
%      abcorr   the scalar string identifying the aberration corrections to
%               apply when computing the observer-target state.  'abcorr'
%               may be any of the following.
%
%                  'NONE'   Apply no correction. Return the
%                           geometric sub-solar point on the target
%                           body.
%
%                  'LT'     Correct for planetary (light time)
%                           aberration.  Both the state and rotation
%                           of the target body are corrected for one
%                           way light time from target to observer.
%
%                           The state of the sun relative to the
%                           target is corrected for one way light
%                           from the sun to the target; this state
%                           is evaluated at the epoch obtained by
%                           retarding et by the one way light time
%                           from target to observer.
%
%                  'LT+S'   Correct for planetary (light time) and
%                           stellar aberrations.  Light time
%                           corrections are the same as in the 'LT'
%                           case above.  The target state is
%                           additionally corrected for stellar
%                           aberration as seen by the observer, and
%                           the sun state is corrected for stellar
%                           aberration as seen from the target.
%
%                  'CN'     Converged Newtonian light time correction.
%                           This option produces a solution that is at
%                           least as accurate at that obtainable 
%                           with the 'LT' option. Whether the 'CN' 
%                           solution is substantially more accurate 
%                           depends on the geometry of the 
%                           participating objects and on the 
%                           accuracy of the input data. In all 
%                           cases this routine will execute more 
%                           slowly when a converged solution is 
%                           computed. See the section titled "The 
%                           Computation of Light Time" in the SPK
%                           Required Reading document spk.req for 
%                           details.
%
%                  'CN+S'   Converged Newtonian light time
%                           and stellar aberration corrections. 
%                           Light time and stellar aberration
%                           corrections are applied as in the
%                           'LT+S' case.
%
%      obsrvr   the scalar sting describing the name of the observing body.
%               This is typically a spacecraft, the earth, or a surface point
%               on the earth.  `obsrvr' is case-insensitive, and leading and
%               trailing blanks in `obsrvr' are not significant.
%               Optionally, you may supply a string containing the
%               integer ID code for the object.  For example both
%               'EARTH' and '399' are legitimate strings that indicate
%               the earth is the observer.
%
%   the call:
%
%      spoint = cspice_subsol( method, target, et, abcorr, obsrvr )
%
%   returns:
%
%      spoint   the double precision 3x1 array or 3xN array describing the
%               sub-solar point on the target body at 'et',
%               expressed relative to the body-fixed frame of the
%               target body.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Find the sub solar position on the earth as seen from the moon at
%      % at epoch JAN 1, 2006 using the 'near point' then the 'intercept'
%      % options. Apply light time correction to return apparent position.
%      %
%      % Load the meta kernel listing the needed SPK, PCK, LSK
%      % kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      et = cspice_str2et( 'JAN 1, 2006' );
%
%      %
%      % First use option 'Near Point'
%      %
%      point1 = cspice_subsol( 'near point', 'earth', et, 'lt+s', 'moon');
%
%      disp( 'Sub solar location coordinates - near point:' )
%      fprintf( '    %15.8f\n', point1 )
%
%      disp(' ')
%
%      %
%      % Now use option 'Intercept'
%      %
%      point2 = cspice_subsol( 'intercept', 'earth', et, 'lt+s', 'moon');
%
%      disp( 'Sub solar location coordinates - intercept' )
%      fprintf( '    %15.8f\n', point2 )
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
%     Sub solar location coordinates - near point:
%          -5872.12723742
%            -91.14115623
%          -2479.72440750
%
%     Sub solar location coordinates - intercept
%          -5866.09275043
%            -91.04749509
%          -2493.87447851
%
%     Distance between locations            (km): 15.38338
%     Angular separation between locations (deg):  0.13826
%
%-Particulars
%
%   cspice_subsol computes the sub-solar point on a target body, as seen by
%   a specified observer.
%
%   There are two different popular ways to define the sub-solar point:
%   "nearest point on target to the sun" or "target surface intercept of
%   line containing target and sun."  These coincide when the target is
%   spherical and generally are distinct otherwise.
%
%   When comparing sub-point computations with results from sources
%   other than SPICE, it's essential to make sure the same geometric
%   definitions are used.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine subsol_c.
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
%      correction.
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
%   -Mice Version 1.0.0, 07-MAR-2007, EDW (JPL)
%
%-Index_Entries
%
%   DEPRECATED sub-solar point
%
%-&

function [spoint] = cspice_subsol( method, target, et, abcorr, obsrvr )

   switch nargin
      case 5

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);

      otherwise

         error ( ['Usage: [_spoint(3)_ ] = '            ...
                  'cspice_subsol( `method`, `target`, ' ...
                                  '_et_, `abcorr`, `obsrvr`)']  )

   end

   %
   % Call the MEX library.
   %
   try
      [spoint] = mice('subsol_c', method, target, et, abcorr, obsrvr);
   catch
      rethrow(lasterror)
   end


