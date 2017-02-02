%-Abstract
%
%   CSPICE_NPELPT calculates the location on an ellipse closest
%   to a specified point, both in three-dimensional space,
%   and the distance between the ellipse and the point.
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
%      point     a 3-vector defining a location in 3-dimensional space.
%
%                [3,1] = size(point); double = class(point)
%
%      ellipse   a SPICE ellipse.
%
%                [1,1] = size(ellipse); struct = class(ellipse)
%
%                The structure has the fields:
%
%                center:    [3,1] = size(center); double = class(center)
%                semiMinor: [3,1] = size(semiMinor); double = class(semiMinor)
%                semiMajor: [3,1] = size(semiMajor); double = class(semiMajor)
%
%   the call:
%
%      [ pnear, dist ] = cspice_npelpt( point, ellipse )
%
%   returns:
%
%      pnear   the 3-vector defining the location on the ellipse nearest
%              to 'point'.
%
%              [3,1] = size(point); double = class(point)
%
%      dist    the distance between the calculated point 'pnear' 
%              and 'point'.
%
%              dist = || pnear - point ||
%
%              [1,1] = size(dist); double = class(dist)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%-Particulars
%
%   Given an ellipse and a point in 3-dimensional space, if the
%   orthogonal projection of the point onto the plane of the ellipse
%   is on or outside of the ellipse, then there is a unique point on
%   the ellipse closest to the original point.  This routine finds
%   that nearest point on the ellipse.  If the projection falls inside
%   the ellipse, there may be multiple points on the ellipse that are
%   at the minimum distance from the original point.  In this case,
%   one such closest point will be returned.
%
%   This routine returns a distance, rather than an altitude, in
%   contrast to the Icy routine cspice_nearpt.  Because our ellipse is
%   situated in 3-space and not 2-space, the input point is not
%   `inside' or `outside' the ellipse, so the notion of altitude does
%   not apply to the problem solved by this routine.  In the case of
%   cspice_nearpt, the input point is on, inside, or outside the ellipsoid,
%   so it makes sense to speak of its altitude.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine npelpt_c.
%
%   MICE.REQ
%   ELLIPSES.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 14-NOV-2014, EDW (JPL)
%
%-Index_Entries
%
%   nearest point on ellipse to point
%
%-&

function [ pnear, dist ] = cspice_npelpt( point, ellipse )

   switch nargin
      case 2

         point  = zzmice_dp(point);
         ellipse= zzmice_ell(ellipse);

      otherwise

         error ( 'Usage: [ pnear, dist ] = cspice_npelpt( point, ellipse )' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [npelpt] = mice( 'npelpt_s', point, ellipse );
      pnear    = reshape( [npelpt.pos], 3, [] );
      dist     = reshape( [npelpt.alt], 1, [] );
   catch
      rethrow(lasterror)
   end


