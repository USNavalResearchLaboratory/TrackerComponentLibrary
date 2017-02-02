%-Abstract
%
%   CSPICE_VPRJP orthogonally projects a vector onto a specified plane.
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
%      vin     the 3-vector to orthogonally project onto a specified plane.
%
%              [3,1] = size(vin); double = class(vin)
%
%      plane   a structure describing a SPICE plane onto which to
%              project 'vin'.
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
%      vout = cspice_vprjp( vin, plane )
%
%   returns:
%
%      vout   3-vector resulting from the orthogonal projection of 'vin'
%             onto 'plane'. 'vout' is the closest point in the specified
%             plane to 'vin'.
%
%             [3,1] = size(vout); double = class(vout)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Find the closest point in the ring plane of a planet to a
%      % spacecraft located at vec1 (in body-fixed coordinates).
%      %
%      vec1 = [ -5., 7., 2.2]';
%
%      %
%      % Define the vector normal as the normal to the
%      % equatorial ring plane, and the origin at body/ring center
%      %
%      norm = [ 0., 0., 1. ]';
%      orig = [ 0., 0., 0. ]';
%
%      %
%      % Create the plane structure.
%      %
%      ring_plane = cspice_nvp2pl( norm, orig );
%
%      %
%      % Project the position vector onto the ring plane...
%      %
%      proj = cspice_vprjp( vec1, ring_plane )
%
%   MATLAB outputs:
%
%      proj =
%
%          -5
%           7
%           0
%
%-Particulars
%
%   Projecting a vector v orthogonally onto a plane can be thought of
%   as finding the closest vector in the plane to v.  This `closest
%   vector' always exists; it may be coincident with the original
%   vector.
%
%   Two related routines are cspice_vprjpi, which inverts an orthogonal
%   projection of a vector onto a plane, and cspice_vproj, which projects
%   a vector orthogonally onto another vector.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine vprjp_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 27-AUG-2012, EDW (JPL)
%
%-Index_Entries
%
%   vector projection onto plane
%
%-&

function [vout] = cspice_vprjp( vin, plane )

   switch nargin

      case 2

         vin   = zzmice_dp( vin );
         plane = zzmice_pln( plane );

      otherwise

         error ( 'Usage: [vout(3)] = cspice_vprjp( vin(3), plane )' )

   end

   %
   % Call the MEX library.
   %
   % The developer decided to not complicate the interface call and so
   % use the individual fields of the 'plane' structure as arguments.
   %
   try
      [vout] = mice('vprjp_c', vin, plane );
   catch
      rethrow(lasterror)
   end

