%-Abstract
%
%   CSPICE_NVC2PL constructs a SPICE plane from a normal vector and a point
%   on the plane.
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
%      normal   [3,1] = size(normal); double = class(normal)
%
%      point    [3,1] = size(point); double = class(point)
%
%               are, respectively, a normal vector and point that
%               define a plane in three-dimensional space.  normal
%               need not be a unit vector. Let the symbol < a, b >
%               indicate the inner product of vectors a and b;
%               then the geometric plane is the set of vectors x
%               in three-dimensional space that satisfy
%
%                  < x - point, normal >  =  0.
%
%   the call:
%
%      plane = cspice_nvp2pl( normal, point )
%
%   returns:
%
%      plane   a structure describing a SPICE plane defined by
%              'normal' and 'point'
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
%      % Define a normal vector from a plane and a
%      % point in a plane.
%      %
%      normal = [ -1.;  5.;   -3.5 ];
%      point  = [  9.; -0.65; -12. ];
%
%      %
%      % Create a plane from the vectors.
%      %
%      plane = cspice_nvp2pl( normal, point );
%
%      %
%      % Calculate a point in the plane, and
%      % two spanning vectors in the plane such that
%      % the point and spanning are mutually orthogonal.
%      %
%      [point, span1, span2] = cspice_pl2psv( plane )
%
%      %
%      % Test point, span1, and span2 orthogonality. The dot
%      % products of any two vectors should equal zero to
%      % within round-off.
%      %
%      fprintf( 'dot( point, span1) = %18.15e\n', dot( point, span1) )
%      fprintf( 'dot( point, span2) = %18.15e\n', dot( point, span2) )
%      fprintf( 'dot( span1, span2) = %18.15e\n', dot( span1, span2) )
%
%    Matlab outputs:
%
%        point =
%
%            -7.777777777777776e-01
%             3.888888888888888e+00
%            -2.722222222222222e+00
%
%
%        span1 =
%
%                                 0
%             5.734623443633283e-01
%             8.192319205190405e-01
%
%
%        span2 =
%
%             9.868415319342446e-01
%             1.324619505952006e-01
%            -9.272336541664042e-02
%
%        dot( point, span1) = 0.000000000000000e+00
%        dot( point, span2) = 5.551115123125783e-17
%        dot( span1, span2) = 0.000000000000000e+00
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
%   the CSPICE routine nvp2pl_c.
%
%   MICE.REQ
%   PLANES.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 27-AUG-2012, EDW (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 30-DEC-2008, EDW (JPL)
%
%-Index_Entries
%
%   normal vector and point to plane
%
%-&

function [plane] = cspice_nvp2pl( normal, point )

   switch nargin

      case 2

         normal = zzmice_dp(normal);
         point  = zzmice_dp(point);

      otherwise

         error ( ['Usage: [plane] = cspice_nvp2pl( normal(3), point(3) )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [plane] = mice('nvp2pl_c', normal, point );
   catch
      rethrow(lasterror)
   end

