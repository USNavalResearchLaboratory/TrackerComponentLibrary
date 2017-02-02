%-Abstract
%
%   CSPICE_PSV2PL returns a SPICE plane given a point and two
%   spanning vectors.
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
%      point    [3,1] = size(point); double = class(point)
%
%      span1    [3,1] = size(span1); double = class(span1)
%
%      span2    [3,1] = size(span2); double = class(span2)
%
%               are, respectively, a point and two spanning vectors
%               that define a geometric plane in three-dimensional
%               space. The plane is the set of vectors
%
%                  point   +   s * span1   +   t * span2
%
%               where 's' and 't' are real numbers.  The spanning
%               vectors 'span1' and 'span2' must be linearly
%               independent, but they need not be orthogonal or
%               unitized.
%
%   the call:
%
%      [plane] = cspice_psv2pl(point, span1, span2 )
%
%   returns:
%
%      plane   a structure describing a SPICE plane defined
%              by 'point', 'span1', and 'span2'.
%
%              [1,1] = size(plane); struct = class(plane)
%
%              The structure has the fields:
%
%                 normal:     [3,1] = size(normal); double = class(normal)
%                 constant:   [1,1] = size(constant); double = class(constant)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Calculate the inclination of the Moon's orbit plane about
%      % the Earth to the orbit plane of the Earth around the sun.
%      %
%      % We want a geometric analysis, so the calculation requires
%      % no aberration correction.
%      %
%      epoch = 'Jan 1 2005';
%      frame = 'ECLIPJ2000';
%      corr  = 'NONE';
%
%      %
%      % Load the kernels we need to retrieve state data.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Convert the time string to  ephemeris time
%      %
%      et = cspice_str2et( epoch );
%
%      %
%      % Calculate the orbit plane of the Earth about
%      % the solar system barycenter at epoch.
%      %
%      [state, ltime] = cspice_spkezr( 'EARTH', et, frame, corr, ...
%                                      'Solar System Barycenter' );
%
%      es_plane            = cspice_psv2pl( state(1:3), ...
%                                           state(1:3), ...
%                                           state(4:6) );
%      [es_norm, es_const] = cspice_pl2nvc( es_plane );
%
%      %
%      % Calculate the orbit plane of the Moon with respect to
%      % the Earth-Moon barycenter at epoch.
%      %
%      [state, ltime] = cspice_spkezr( 'MOON', et, frame, corr, ...
%                                      'EARTH BARYCENTER' );
%
%      em_plane            = cspice_psv2pl( state(1:3), ...
%                                           state(1:3), ...
%                                           state(4:6) );
%      [em_norm, em_const] = cspice_pl2nvc( em_plane );
%
%      %
%      % Calculate the inclination equals (output in degrees).
%      % Depending on the orientation of the plane normals, the
%      % cspice_vsep result may exceed 90 degrees. If, so subtract
%      % the value off 180 degrees.
%      %
%      loc_inc = cspice_vsep( es_norm, em_norm );
%
%      if ( loc_inc > cspice_halfpi )
%         loc_inc = cspice_pi - loc_inc;
%      end
%
%      fprintf( 'Moon-Earth orbit plane inclination (degrees): %f\n', ...
%             loc_inc * cspice_dpr )
%
%   Matlab outputs:
%
%      Moon-Earth orbit plane inclination (degrees): 5.042496
%
%-Particulars
%
%   Mice geometry routines that deal with planes use the `plane'
%   data type to represent input and output planes.  This data type
%   makes the subroutine interfaces simpler and more uniform.
%
%   The Mice routines that produce SPICE planes from data that
%   define a plane are:
%
%      cspice_nvc2pl ( Normal vector and constant to plane )
%      cspice_nvp2pl ( Normal vector and point to plane    )
%      cspice_psv2pl ( Point and spanning vectors to plane )
%
%   The Mice routines that convert SPICE planes to data that
%   define a plane are:
%
%      cspice_pl2nvc ( Plane to normal vector and constant )
%      cspice_pl2nvp ( Plane to normal vector and point    )
%      cspice_pl2psv ( Plane to point and spanning vectors )
%
%   Any of these last three routines may be used to convert this
%   routine's output, 'plane', to another representation of a
%   geometric plane.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine psv2pl_c.
%
%   MICE.REQ
%   PLANES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 27-AUG-2012, EDW (JPL)
%
%-Index_Entries
%
%   plane to point and spanning vectors
%
%-&

function [plane] = cspice_psv2pl( point, span1, span2 )

   switch nargin

      case 3

         point = zzmice_dp(point);
         span1 = zzmice_dp(span1);
         span2 = zzmice_dp(span2);

      otherwise

         error ( ['Usage: [plane] = ' ...
                  'cspice_psv2pl( point(3), span1(3), span2(3) )'] )

   end

   %
   % Call the MEX library.
   %
   % The developer decided to not complicate the interface call and so
   % use the individual fields of the 'plane' structure as arguments.
   %
   try
      [plane] = mice('psv2pl_c', point, span1, span2 );
   catch
      rethrow(lasterror)
   end

