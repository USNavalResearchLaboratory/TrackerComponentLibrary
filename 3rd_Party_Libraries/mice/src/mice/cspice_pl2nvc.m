%-Abstract
%
%   CSPICE_PL2NVC returns a unit normal vector and constant defining
%   a specified plane.
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
%      plane   a structure describing a SPICE plane.
%
%              [1,1] = size(plane); struct = class(plane)
%
%              The structure has the fields:
%
%                 normal:     [3,1] = size(normal); double = class(normal)
%                 constant:   [1,1] = size(constant); double = class(constant)
%
%   the call:
%
%      [normal, constant] = cspice_pl2nvc( plane )
%
%   returns:
%
%      normal     [3,1] = size(normal); double = class(normal)
%
%      constant   [1,1] = size(constant); double = class(constant)
%
%                 are, respectively, a unit normal vector and
%                 constant that define the geometric plane
%                 represented by 'plane'.  Let the symbol < a, b >
%                 indicate the inner product of vectors a and b; then
%                 the geometric plane is the set of vectors x in
%                 three-dimensional space that satisfy
%
%                    < x,  normal >  =  constant.
%
%                 'normal' is a unit vector. 'constant' is the distance
%                 of the plane from the origin;
%
%                    constant * normal
%
%                 is the closest point in the plane to the origin.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % A simple task, determine the distance of a plane
%      % from the origin.
%      %
%      % Define the plane with a vector normal to the plane
%      % and a point in the plane.
%      %
%      normal = [ -1.;  5.;    -3.5 ];
%      point  = [  9.; -0.65;  -12. ];
%
%      %
%      % create the SPICE plane from the normal and constant.
%      %
%      plane = cspice_nvp2pl( normal, point );
%
%      %
%      % Calculate the normal vector and constant defining
%      % the plane. The constant value is the distance from
%      % the origin to the plane.
%      %
%      [normal, constant ] = cspice_pl2nvc( plane )
%
%      %
%      % Confirm the results. Calculate a vector
%      % from the origin to the plane.
%      %
%      vec = constant * normal;
%
%      %
%      % Now calculate a vector in the plane from the
%      % location in the plane defined by 'vec'.
%      %
%      plane_vec = vec - point;
%
%      %
%      % These vectors should be orthogonal.
%      %
%      dot( plane_vec, vec )
%
%   MATLAB outputs:
%
%      normal =
%
%          -1.616904166908886e-01
%           8.084520834544432e-01
%          -5.659164584181102e-01
%
%
%      constant =
%
%           4.810289896553937e+00
%
%   The dot product result to check orthogonality...
%
%      ans =
%
%          -3.552713678800501e-15
%
%   Zero, to double precision round-off, so orthogonal to that
%   precision.
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
%   the CSPICE routine pl2nvc_c.
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
%   plane to normal vector and constant
%
%-&

function [normal, constant] = cspice_pl2nvc( plane )

   switch nargin

      case 1

         plane = zzmice_pln( plane );

      otherwise

         error( ['Usage: [normal(3), constant] = cspice_pl2nvc( plane )'] )

   end

   %
   % Call the MEX library.
   %
   % The developer decided to not complicate the interface call and so
   % use the individual fields of the 'plane' structure as arguments.
   %
   try
      [normal, constant] = mice('pl2nvc_c', plane );
   catch
      rethrow(lasterror)
   end



