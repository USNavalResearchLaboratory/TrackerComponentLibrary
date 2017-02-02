%-Abstract
%
%   CSPICE_PL2NVP returns a unit normal vector and point that define
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
%      [normal, point] = cspice_pl2nvp( plane )
%
%   returns:
%
%      normal   [3,1] = size(normal); double = class(normal)
%
%      point    [3,1] = size(point); double = class(point)
%
%               are, respectively, a unit normal vector and point
%               that define the geometric plane represented by
%               plane.  Let the symbol < a, b > indicate the inner
%               product of vectors a and b; then the geometric
%               plane is the set of vectors x in three-dimensional
%               space that satisfy
%
%                  < x - point, normal >  =  0.
%
%               'point' is always the closest point in the input
%               plane to the origin. 'point' is always a
%               non-negative scalar multiple of 'normal'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % A trivial example of plane creation and
%      % decomposition. Create a plane via the definition
%      % of a plane normal and constant, determine a
%      % point in the plane.
%      %
%      plane_norm = [ 2.44; -5./3.; 11./9. ];
%      const      = cspice_pi;
%
%      %
%      % Construct the SPICE plane.
%      %
%      plane = cspice_nvc2pl( plane_norm, const );
%
%      %
%      % Convert the plane to a normal vector, point
%      % representation, 'point' lies in the plane.
%      %
%      [ norm_vec, point ] = cspice_pl2nvp( plane )
%
%   MATLAB outputs:
%
%      norm_vec =
%
%           7.630514391556131e-01
%          -5.212099994232330e-01
%           3.822206662437042e-01
%
%
%      point =
%
%           7.496657642594705e-01
%          -5.120667788657586e-01
%           3.755156378348896e-01
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
%   the CSPICE routine pl2nvp_c.
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
%   plane to normal vector and point
%
%-&

function [normal, point] = cspice_pl2nvp( plane )

   switch nargin

      case 1

         plane = zzmice_pln( plane );

      otherwise

         error ( ['Usage: [normal(3), point(3)] = cspice_pl2nvp( plane )'] )

   end

   %
   % Call the MEX library.
   %
   % The developer decided to not complicate the interface call and so
   % use the individual fields of the 'plane' structure as arguments.
   %
   try
      [normal, point] = mice('pl2nvp_c', plane );
   catch
      rethrow(lasterror)
   end

