%-Abstract
%
%   CSPICE_EUL2M constructs a 3x3, double precision rotation matrix
%   from a set of Euler angles and the corresponding rotation axes.
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
%      angle3
%      angle2
%      angle1    a set of scalar or 1xN arrays of double precision
%                rotation angles measured in radians
%      axis3
%      axis2
%      axis1     the scaler integer indices defining the rotation axis
%                corresponding to each angle
%
%                The values of axisX may be 1, 2, or 3, indicating
%                the x, y, and z axes respectively.
%
%   the call:
%
%      r = cspice_eul2m(angle3, angle2, angle1, axis3, axis2, axis1)
%
%   returns:
%
%      r   a 3x3 or 3x3xN array of double precision matrices defined by
%          the Euler rotation scheme, with r defined as:
%
%          r = [ angle3 ]     [ angle2 ]      [ angle1 ]
%                        axis3          axis2           axis1
%
%         'r' return with the same vectorization measure (N) as
%         'angle3', 'angle2', and 'angle1'.
%
%      Note: the rotation defines a coordinate system rotation,
%      e.g. a single rotation of 90 degrees about Z maps
%      the vector [ 1, 0, 0] (the +x unit vector) to [ 0, -1, 0]
%      (the -y unit vector). A vector rotation would map the
%      +x unit vector to the +y unit vector.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Create the rotation matrix for a single coordinate
%      % rotation of 90 degrees about the Z axis. As the
%      % second and third angles are 0, the final two axes IDs,
%      % 1, 1, have no effect for in this example.
%      %
%      rot = cspice_eul2m( cspice_halfpi, 0, 0, 3, 1, 1 );
%
%      %
%      % Output the result of rotating the +x unit vector
%      % using the 'rot' matrix.
%      %
%      rot * [1; 0; 0 ]
%
%   MATLAB outputs:
%
%      ans =
%
%          0.0000
%         -1.0000
%               0
%
%   A representation of the -y unit vector to round-off error accuracy.
%
%-Particulars
%
%   A word about notation:  the symbol
%
%      [ x ]
%           i
%
%   indicates a rotation of x radians about the ith coordinate axis.
%   To be specific, the symbol
%
%      [ x ]
%           1
%
%   indicates a coordinate system rotation of x radians about the
%   first, or x-, axis; the corresponding matrix is
%
%      +-                    -+
%      |  1      0       0    |
%      |                      |
%      |  0    cos(x)  sin(x) |.
%      |                      |
%      |  0   -sin(x)  cos(x) |
%      +-                    -+
%
%   Remember, this is a COORDINATE SYSTEM rotation by x radians; this
%   matrix, when applied to a vector, rotates the vector by -x
%   radians, not x radians.  Applying the matrix to a vector yields
%   the vector's representation relative to the rotated coordinate
%   system.
%
%   The analogous rotation about the second, or y-, axis is
%   represented by
%
%      [ x ]
%           2
%
%   which symbolizes the matrix
%
%      +-                    -+
%      | cos(x)   0   -sin(x) |
%      |                      |
%      |  0       1      0    |,
%      |                      |
%      | sin(x)   0    cos(x) |
%      +-                    -+
%
%   and the analogous rotation about the third, or z-, axis is
%   represented by
%
%      [ x ]
%           3
%
%   which symbolizes the matrix
%
%      +-                    -+
%      |  cos(x)  sin(x)   0  |
%      |                      |
%      | -sin(x)  cos(x)   0  |.
%      |                      |
%      |  0        0       1  |
%      +-                    -+
%
%   From time to time, (depending on one's line of work, perhaps) one
%   may happen upon a pair of coordinate systems related by a
%   sequence of rotations.  For example, the coordinate system
%   defined by an instrument such as a camera is sometime specified
%   by RA, DEC, and twist; if alpha, delta, and phi are the rotation
%   angles, then the series of rotations
%
%      [ phi ]     [ pi/2  -  delta ]     [ alpha ]
%             3                      2             3
%
%   produces a transformation from inertial to camera coordinates.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine eul2m_c.
%
%   MICE.REQ
%   ROTATION.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   euler angles to matrix
%
%-&

function [r] = cspice_eul2m(angle3, angle2, angle1, axis3, axis2, axis1)

   switch nargin
      case 6

         angle3 = zzmice_dp(angle3);
         angle2 = zzmice_dp(angle2);
         angle1 = zzmice_dp(angle1);
         axis3  = zzmice_int(axis3);
         axis2  = zzmice_int(axis2);
         axis1  = zzmice_int(axis1);

      otherwise

         error( ['Usage: [_r(3,3)_] = '                         ...
                 'cspice_eul2m(_angle3_, _angle2_, _angle1_, '  ...
                 'axis3, axis2, axis1)']  )

   end

   %
   % Call the MEX library.
   %
   try
      [r] = mice('eul2m_c',angle3,angle2,angle1,axis3,axis2,axis1);
   catch
      rethrow(lasterror)
   end




