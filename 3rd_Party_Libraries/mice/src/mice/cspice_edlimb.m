%-Abstract
%
%   CSPICE_EDLIMB calculates the limb of a triaxial ellipsoid
%   as viewed from a specified location.
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
%      a
%      b
%      c        [1,1] = size(a); double = class(a)
%               [1,1] = size(b); double = class(b)
%               [1,1] = size(c); double = class(c)
%
%               are the  lengths of the semi-axes of a triaxial ellipsoid.
%               The ellipsoid is centered at the origin and oriented so that
%               its axes lie on the x, y and z axes. 'a', 'b', and 'c' are
%               the lengths of the semi-axes that respectively point in the
%               x, y, and z directions.
%
%      viewpt   a point from which the ellipsoid is viewed. 'viewpt' must be
%               outside of the ellipsoid.
%
%               [3,1] = size(viewpt); double = class(viewpt)
%
%   the call:
%
%      limb = cspice_edlimb( a, b, c, viewpt )
%
%   returns:
%
%      limb   the SPICE ellipse that represents the limb of the ellipsoid
%             observed from 'viewpt'.
%
%              [1,1] = size(limb); struct = class(limb)
%
%              The structure has the fields:
%
%                 center:    [3x1 double]
%                 semiMajor: [3x1 double]
%                 semiMinor: [3x1 double]
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Define an ellipsoid
%      %
%      a = sqrt(2.);
%      b = 2.*sqrt(2.);
%      c = sqrt(2.);
%
%      %
%      % Locate a viewpoint exterior to the ellipsoid.
%      %
%      viewpt = [ 2., 0.,  0. ]';
%
%      %
%      % Calculate the limb ellipse as seen by from the viewpoint.
%      %
%      limb = cspice_edlimb( a, b, c, viewpt );
%
%      %
%      % Output the structure components.
%      %
%      smin   = limb.semiMinor
%      smaj   = limb.semiMajor
%      center = limb.center
%
%      %
%      % Check against expected values:
%      %
%      % Semiminor: 0., 0., -1.
%      % Semimajor: 0., 2.,  0.
%      % Center   : 1., 0.,  0.
%      %
%
%   MATLAB outputs:
%
%      smin =
%
%           0
%           0
%          -1
%
%
%      smaj =
%
%           0
%           2
%           0
%
%
%      center =
%
%           1.000000000000000e+00
%                               0
%                               0
%
%-Particulars
%
%   The limb of a body, as seen from a viewing point, is the boundary
%   of the portion of the body's surface that is visible from that
%   viewing point.  In this definition, we consider a surface point
%   to be `visible' if it can be connected to the viewing point by a
%   line segment that doesn't pass through the body.  This is a purely
%   geometrical definition that ignores the matter of which portions
%   of the surface are illuminated, or whether the view is obscured by
%   any additional objects.
%
%   If a body is modeled as a triaxial ellipsoid, the limb is always
%   an ellipse.  The limb is determined by its center, a semi-major
%   axis vector, and a semi-minor axis vector.
%
%   We note that the problem of finding the limb of a triaxial
%   ellipsoid is mathematically identical to that of finding its
%   terminator, if one makes the simplifying assumption that the
%   terminator is the limb of the body as seen from the vertex of the
%   umbra.  So, this routine can be used to solve this simplified
%   version of the problem of finding the terminator.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine edlimb_c.
%
%   MICE.REQ
%   ELLIPSES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 09-NOV-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   ellipsoid limb
%
%-&

function [ limb ] = cspice_edlimb( a, b, c, viewpt )

   switch nargin
      case 4

         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);
         viewpt = zzmice_dp(viewpt);

      otherwise

         error ( 'Usage: [ limb ] = cspice_edlimb( a, b, c, viewpt )' )

   end

   try
      [ limb ] = mice( 'edlimb_s', a, b, c, viewpt );
   catch
      rethrow(lasterror)
   end


