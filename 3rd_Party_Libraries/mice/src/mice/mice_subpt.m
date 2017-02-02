%-Abstract
%
%   MICE_SUBPT determines the coordinates of the sub-observer point
%   on a target body at a particular epoch, optionally corrected
%   for planetary (light time) and stellar aberration. The call also
%   returns the observer's altitude above the target body.
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
%      method   a string defining the method to use in the calculation
%
%      target   the string name of the observed target body
%
%      et       the double precision ephemeris time of interest
%
%      abcorr   the string defining the aberration correction to use in
%               the calculation
%
%      obsrvr   the string name of the observing body
%
%   the call:
%
%      [spoint] = mice_subpt( method, target, et, abcorr, obsrvr)
%
%   returns:
%
%      spoint   the scalar or 1xN array of structures, each structure
%               consisting of two fields:
%
%                  'pos'   the double-precision 3-vector containing the
%                          coordinates of the 'obsrvr' subpoint on 'target'
%                          relative to the body-fixed frame of 'target'
%
%                  'alt'   the double precision scalar altitude of 'obsrvr'
%                          above 'target'
%
%              'spoint' returns with the same vectorization measure (N)
%               as 'et'.
%
%      Note, If needed the user can extract the field data from vectorized
%      'spoint' structures into separate arrays.
%
%      Extract the N 'pos' field data to a 3XN array 'position':
%
%         position = reshape( [spoint(:).pos], 3, [] )
%
%      Extract the N 'alt' field data to a 1XN array 'altitude':
%
%         altitude = reshape( [point(:).alt], 1, [] )
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
%      % a given time using the "near point" then the "intercept"
%      % options.
%      %
%      % Load the meta kernel listing the needed SPK, PCK, LSK
%      % kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Calculate the location of the sub point of the moon as
%      % seen from Earth at epoch JAN 1, 2006. Apply light time
%      % correction to return apparent position.
%      %
%      et = cspice_str2et( 'JAN 1, 2006' );
%
%      %
%      % First use option 'Near Point'
%      %
%      [point1] = mice_subpt( 'near point', 'earth', et, 'lt+s', 'moon');
%
%      %
%      % Now use option 'Intercept'
%      %
%      [point2] = mice_subpt( 'intercept', 'earth', et, 'lt+s', 'moon');
%
%      %
%      % Calculate the Euclidean distance between the two locations
%      % and the angular separation between the position vectors.
%      %
%      dist = norm( point1.pos - point2.pos);
%      sep  = cspice_vsep(point1.pos, point2.pos )*cspice_dpr;
%
%      txt = sprintf( 'Distance between locations            (km): %8.5f', ...
%                                                                       dist);
%      disp( txt )
%
%      txt = sprintf( 'Angular separation between locations (deg): %8.5f', ...
%                                                                       sep );
%      disp( txt )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Distance between locations            (km): 16.70963
%      Angular separation between locations (deg):  0.15020
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
%      [point] = mice_subpt( 'near point', 'earth', et, 'lt+s', 'moon');
%
%      %
%      % Convert the subpoint coordinates to lat/lon expressed in degrees with
%      % the radius.
%      %
%      % Extract from the 'point' structure the 3XN array of position data.
%      %
%      position = reshape( [point(:).pos], 3, [] )
%
%      [radius, longitude, latitude] = cspice_reclat(position);
%      longitude                     = longitude * cspice_dpr;
%      latitude                      = latitude  * cspice_dpr;
%
%      %
%      % Convert the 'et' epochs to calendar format.
%      %
%      utc = cspice_et2utc( et, 'C', 3 );
%
%      for i=1:13
%         txt = sprintf( 'Moon subpoint epoch: %s', utc(i,:) );
%         disp( txt )
%
%         txt = sprintf( '              (deg): longitude %8.4f', longitude(i) );
%         disp( txt )
%
%         txt = sprintf( '              (deg): latitude  %8.4f', latitude(i) );
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
%   A sister version of this routine exists named cspice_subpt that returns
%   the structure field data as separate arguments.
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
%   -Mice Version 1.0.0, 16-DEC-2005, EDW (JPL)
%
%-Index_Entries
%
%   sub-observer point
%
%-&

function [spoint] = mice_subpt( method, target, et, abcorr, obsrvr )

   switch nargin
      case 5

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);

      otherwise
         error ( ['Usage: [_spoint_] = '                            ...
                  'mice_subpt( `method`, `target`, _et_, `abcorr`, `obsrvr`)'])
   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [spoint] = mice('subpt_s', method, target, et, abcorr, obsrvr);
   catch
      rethrow(lasterror)
   end


