%-Abstract
%
%   CSPICE_INELPL finds the intersection of an ellipse and a plane.
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
%      ellipse   a structure describing a SPICE ellipse. The ellipse is
%                allowed to be degenerate: one or both semi-axes may
%                have zero length.
%
%                [1,1] = size(ellipse); struct = class(ellipse)
%
%                The structure has the fields:
%
%                  center:    [3,1] = size(center); double = class(center)
%                  semiMajor: [3,1] = size(semiMajor); double = class(semiMajor)
%                  semiMinor: [3,1] = size(semiMinor); double = class(semiMinor)
%
%      plane     a structure describing a SPICE plane. The intersection of
%                'plane' and 'ellipse' is sought.
%
%                [1,1] = size(plane); struct = class(plane)
%
%                The structure has the fields:
%
%                   normal:     [3,1] = size(normal); double = class(normal)
%                   constant:   [1,1] = size(constant); double = class(constant)
%
%   the call:
%
%      [ nxpts, xpt1, xpt2] = cspice_inelpl( ellipse, plane )
%
%   returns:
%
%      nxpts   the number of points of intersection of the geometric plane
%              and ellipse represented by `plane' and `ellips'. `nxpts' may
%              take the values 0, 1, 2 or -1. The value -1 indicates the
%              ellipse lies in the plane, so the number of intersection
%              points is infinite.
%
%              -1 also signals for the degenerate case where the ellipse
%              structure defines a single point and that point lies
%              in the plane of interest. In this case, -1 means not an
%              infinite number of intersections, rather that the
%              ellipse is a subset of the plane, that subset having
%              measure one.
%
%              [1,1] = size(nxpts); int32 = class(nxpts)
%
%      xpt1,
%      xpt2    the points of intersection of the input plane and ellipse.
%              If there is only one intersection point, both 'xpt1' and
%              'xpt2' contain that point. If the number of intersection
%              points is zero or infinite, the contents of 'xpt1' and
%              'xpt2' are undefined.
%
%              [3,1] = size(xpt1); double = class(xpt1)
%              [3,1] = size(xpt2); double = class(xpt2)
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
%            KERNELS_TO_LOAD = ( '/kernels/gen/lsk/naif0009.tls'
%                                '/kernels/gen/spk/de421.bsp'
%                                '/kernels/gen/pck/pck00009.tpc'
%                      )
%
%         \begintext
%
%
%      %
%      % Standard SPK, LSK, PCK files.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Retrieve the triaxial radii of Saturn (699)
%      %
%      radii = cspice_bodvrd( 'SATURN', 'RADII', 3 );
%
%      %
%      % Define a position in the body frame at one hundred equatorial
%      % radii out along the x axis, one hundred radii above the
%      % equatorial plane.
%      %
%      vertex = [ 100.0 * radii(1), 0.0, radii(1) *100.0 ]';
%
%      %
%      % Find the limb of the ellipsoid as seen from the
%      % point 'vertex'. 'limb' returns as a CSPICE_ELLIPSE.
%      %
%      limb = cspice_edlimb( radii(1), radii(2), radii(3), vertex );
%
%      %
%      % Define the equatorial plane as a SPICE plane. The Z
%      % axis is normal to the plane, the origin lies in the
%      % plane.
%      %
%      normal = [ 0., 0., 1.]';
%      point  = [ 0., 0., 0.]';
%      plane  = cspice_nvp2pl( normal, point);
%
%      %
%      % Calculate the intersection of the 'limb' and 'plane'.
%      %
%      [ nxpts, xpt1, xpt2] = cspice_inelpl( limb, plane );
%
%      fprintf(                                                             ...
%         'Observer at (100, 0, 100) radii, no. intersection points: %d\n', ...
%          nxpts )
%      fprintf( '   Intersection points\n' )
%      fprintf( '%.12g  %.12g  %.12g\n',   xpt1   )
%      fprintf( '%.12g  %.12g  %.12g\n\n', xpt2   )
%
%      %
%      % One hundred radii along the Z pole vector (positive)
%      %
%      vertex = [ 0.0 * radii(1), 0.0, radii(1) *100.0 ]';
%
%      %
%      % The resulting limb ellipse should lie parallel to, but
%      % not in the same plane as the equatorial plane. No
%      % intersection should exist.
%      %
%      limb = cspice_edlimb( radii(1), radii(2), radii(3), vertex );
%      [ nxpts, xpt1, xpt2] = cspice_inelpl( limb, plane );
%
%      fprintf(                                                            ...
%         'Ellipse/plane parallel case, no. intersection points: %d\n\n',  ...
%         nxpts )
%
%      %
%      % One radii along the X axis, i.e. on the surface, a very
%      % degenerate case.
%      %
%      vertex = [ radii(1), 0.0, 0.0 ]';
%
%      %
%      % In this case the limb ellipse exists as a point at (x, 0, 0).
%      %
%      limb = cspice_edlimb( radii(1), radii(2), radii(3), vertex );
%
%      %
%      % Calculate the intersection of the plane and the degenerate ellipse.
%      %
%      [ nxpts, xpt1, xpt2 ] = cspice_inelpl( limb, plane );
%
%      %
%      % As the point (x, 0, 0) exists in 'plane' and that point represents
%      % the complete ellipse, the routine should return -1 for infinite
%      % number of intersections - though in this case the intersection
%      % contains only one element.
%      %
%      fprintf( 'Degenerate case, no. intersection points: %d\n', nxpts )
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in IDL due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs for the first example::
%
%      Observer at (100, 0, 100) radii, no. intersection points: 2
%         Intersection points
%      602.68   60264.9865247   3.18323145621e-12
%      602.68  -60264.9865247  -9.37916411203e-12
%
%   MATLAB outputs for the second example, we expect no intersection:
%
%      Ellipse/plane parallel case, no. intersection points: 0
%
%   MATLAB outputs for the degenerate case:
%
%      Degenerate case, no. intersection points: 1
%
%-Particulars
%
%   This routine computes the intersection set of a non-degenerate
%   plane with a possibly degenerate ellipse. The ellipse is allowed
%   to consist of a line segment or a point.
%
%   A plane may intersect an ellipse in 0, 1, 2, or infinitely many
%   points. For there to be an infinite set of intersection points,
%   the ellipse must lie in the plane and consist of more than one
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine inelpl_c.
%
%   MICE.REQ
%   ELLIPSES.REQ
%   PLANES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 27-AUG-2012, EDW (JPL)
%
%-Index_Entries
%
%   intersection of ellipse and plane
%
%-&

function [ nxpts, xpt1, xpt2] = cspice_inelpl( ellipse, plane )

   switch nargin
      case 2

         ellipse = zzmice_ell(ellipse);
         plane   = zzmice_pln(plane);

      otherwise

         error ( ['Usage: [ nxpts, xpt1, xpt2] = ' ...
                         ' cspice_inelpl( ellipse, plane )'] )

   end

   try
      [ nxpts, xpt1, xpt2] = mice( 'inelpl_c', ellipse, plane );
   catch
      rethrow(lasterror)
   end


