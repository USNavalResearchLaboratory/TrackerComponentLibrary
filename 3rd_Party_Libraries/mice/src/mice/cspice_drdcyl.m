%-Abstract
%
%   CSPICE_DRDCYL computes the Jacobian of the transformation from
%   cylindrical to rectangular coordinates.
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
%      r     distance of the point of interest from z axis.
%
%            [1,n] = size(r); double = class(r)
%
%      lon   cylindrical angle (in radians) of the point of interest from the xz
%            plane. The angle increases in the counterclockwise sense about the
%            +z axis.
%
%            [1,n] = size(lon); double = class(lon)
%
%      z     height of the point above xy plane.
%
%            [1,n] = size(z); double = class(z)
%
%   the call:
%
%      jacobi = cspice_drdcyl( r, lon, z)
%
%   returns:
%
%      jacobi   the matrix of partial derivatives of the conversion between
%               cylindrical and rectangular coordinates. It has the form
%
%               If [1,1] = size(r) then [3,3]   = size(jacobi)
%               If [1,n] = size(r) then [3,3,n] = size(jacobi)
%                                        double = class(jacobi)
%
%                   -                               -
%                  |  dx/dr     dx/dlon     dx/dz    |
%                  |                                 |
%                  |  dy/dr     dy/dlon     dy/dz    |
%                  |                                 |
%                  |  dz/dr     dz/dlon     dz/dz    |
%                   -                               -
%
%               evaluated at the input values of r, lon and z.  Here x,y, and
%               z are given by the familiar formulae
%
%                  x = r*cos(lon)
%                  y = r*sin(lon)
%                  z = z
%
%-Examples
%
%   None.
%
%-Particulars
%
%   It is often convenient to describe the motion of an object in
%   the cylindrical coordinate system.  However, when performing
%   vector computations its hard to beat rectangular coordinates.
%
%   To transform states given with respect to cylindrical coordinates
%   to states with respect to rectangular coordinates, one uses
%   the Jacobian of the transformation between the two systems.
%
%   Given a state in cylindrical coordinates
%
%      ( r, lon, z, dr, dlon, dz )
%
%   the velocity in rectangular coordinates is given by the matrix
%   equation:
%                  t          |                          t
%      (dx, dy, dz)   = jacobi|          * (dr, dlon, dz)
%                             |(r,lon,z)
%
%   This routine computes the matrix
%
%            |
%      jacobi|
%            |(r,lon,z)
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine drdcyl_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 09-NOV-2012, EDW (JPL), SCK (JPL)
%
%-Index_Entries
%
%   Jacobian of rectangular w.r.t. cylindrical coordinates
%
%-&

function [jacobi] = cspice_drdcyl( r, lon, z)

   switch nargin
      case 3

         r   = zzmice_dp(r);
         lon = zzmice_dp(lon);
         z   = zzmice_dp(z);

      otherwise

         error( 'Usage: [_jacobi(3,3)_] = cspice_drdcyl( _r_, _lon_, _z_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [jacobi] = mice('drdcyl_c', r, lon, z);
   catch
      rethrow(lasterror)
   end




