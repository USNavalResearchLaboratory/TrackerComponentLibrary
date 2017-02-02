%-Abstract
%
%   CSPICE_M2EUL factors a rotation matrix into a product of
%   three rotations about specified coordinate axes.
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
%      r       is a 3x3 or 3x3xN array of rotation matrices to factor as
%              a product of three rotations about a specified
%              coordinate axes.  The angles of these rotations are
%              called "Euler angles".
%
%      axis3
%      axis2
%      axis1   are the scalar integer indices of the rotation axes of the
%              "factor" rotations, whose product is 'r'. 'r' is
%              factored as
%
%                 r = [ angle3 ]     [ angle2 ]     [ angle1 ]
%                               axis3          axis2          axis1
%
%              The axis numbers must belong to the set {1, 2, 3}.
%              The second axis number MUST differ from the first
%              and third axis numbers.
%
%   the call:
%
%      [angle3, angle2, angle1] = cspice_m2eul( r, axis3, axis2, axis1 )
%
%   returns:
%
%      angle3
%      angle2
%      angle1   the scalar or 1xN arrays of double precision
%               Euler angles measured where the angle satisfy
%
%                   r = [ angle3 ]     [ angle2 ]     [ angle1 ]
%                                axis3           axis2          axis1
%
%                  The range of 'angle3' and 'angle1' is (-pi, pi].
%
%                  The range of 'angle2' depends on the exact set of
%                  axes used for the factorization.  For
%                  factorizations in which the first and third axes
%                  are the same,
%
%                     r = [R]  [S]  [T]
%                            a    b    a
%
%                  the range of 'angle2' is [0, pi].
%
%                  For factorizations in which the first and third
%                  axes are different,
%
%                     r = [R]  [S]  [T] ,
%                            a    b    c
%
%                  the range of angle2 is [-pi/2, pi/2].
%
%                  For rotations such that 'angle3' and 'angle1' are not
%                  uniquely determined, 'angle3' will always be set to
%                  zero; 'angle1' is then uniquely determined.
%
%               'angle3', 'angle2', and 'angle1' return with the same
%               vectorization measure (N) as 'r'.
%
%      Note, the call sequence:
%
%         [angle3, angle2, angle1] = cspice_m2eul(r, axis3, axis2, axis1)
%         r = cspice_eul2m(angle3, angle2, angle1, axis3, axis2, axis1)
%
%      preserves 'r' to round-off error.
%
%      Yet, the call sequence:
%
%         r = cspice_eul2m(angle3, angle2, angle1, axis3, axis2, axis1)
%        [angle3, angle2, angle1] = cspice_m2eul(r, axis3, axis2, axis1)
%
%      preserves 'angle3', 'angle2', and 'angle1' only if the initial
%      values of the angle existed within the range of cspice_m2eul's
%      output.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Example(1):
%
%      %
%      % Scalar example, conversion of instrument pointing from a matrix
%      % representation to Euler angles:
%      %
%      % Suppose we want to find camera pointing in 'alpha', 'delta', and
%      % 'kappa', given the inertial-to-camera coordinate transformation
%      %
%      ticam = [                                                            ...
%         [ 0.49127379678135830  0.50872620321864170  0.70699908539882417 ]
%         [ -0.50872620321864193 -0.49127379678135802  0.70699908539882428]
%         [  0.70699908539882406 -0.70699908539882439  0.01745240643728360] ];
%
%      %
%      % We want to find angles alpha, delta, kappa such that
%      %
%      %
%      %   ticam  =  [ kappa ]  [ pi/2 - delta ]  [ pi/2 + alpha ] .
%      %                      3                 1                 3
%      %
%      %
%      % Factor the matrix to the Euler angles corresponding to a
%      % 3-1-3 rotation.
%      %
%      [ kappa, ang2, ang1  ] = cspice_m2eul( ticam, 3, 1, 3 );
%
%      alpha  =  ang1          - cspice_halfpi;
%      delta  =  cspice_halfpi - ang2;
%
%      %
%      %  calculates the desired angles.  If we wish to make sure that
%      % alpha, delta, and kappa are in the ranges [0, 2pi),
%      % [-pi/2, pi/2], and [0, 2pi) respectively, we may add the code
%      %
%
%      if ( alpha < 0. )
%       alpha = alpha + cspice_twopi;
%      end
%
%      if ( kappa < 0. )
%         kappa = kappa + cspice_twopi;
%      end
%
%      %
%      % Output the 3-1-3 Euler rotation angles corresponding to 'ticam'.
%      %
%      fprintf( '%12.5f   %12.5f   %12.5f\n', ...
%               [ alpha, delta, kappa ] *cspice_dpr)
%
%   MATLAB outputs:
%
%         315.00000        1.00000       45.00000
%
%      Example(2):
%
%      %
%      % Vectorized example, input an array of ephemeris times, calculate
%      % the corresponding J2000 to IAU_MOON transformation matrices.
%      %
%      cspice_furnsh('standard.tm')
%
%      et0 = cspice_str2et( 'Jan 1 2000 12:00:00 TDB' );
%      et1 = cspice_str2et( 'Jan 1 2010 12:00:00 TDB' );
%
%      n     = 10;
%      times = et0 + (1:n)* (et1 - et0)/n;
%      quot   = cspice_pxform( 'J2000', 'IAU_MOON', times );
%
%      %
%      % Factor the matrices to the Euler angles corresponding to a
%      % 3-2-1 rotation set.
%      %
%      [a3,a2,a1] = cspice_m2eul( quot, 1,2,3);
%
%      %
%      % Output the 3-2-1 Euler rotation angles corresponding to 'quot'.
%      %
%      fprintf( '%12.5f   %12.5f   %12.5f\n', [a1; a2; a3] * cspice_dpr )
%
%      cspice_kclear
%
%   MATLAB outputs:
%
%         -52.93007       18.11962       15.07397
%          77.30266      -22.59555        3.51974
%        -150.68645       12.42680      -18.79120
%         -14.28248        4.91714       21.55874
%         120.06957      -19.09792      -11.00536
%        -109.73801       20.66329       -7.52692
%          23.54335       -8.43440       20.49467
%         160.13917       -9.11890      -20.58629
%         -66.71201       21.70068        7.52880
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine m2eul_c.
%
%   MICE.REQ
%   ROTATION.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 30-DEC-2008, EDW (JPL)
%
%      Corrected misspellings.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   matrix to euler angles
%
%-&

function [angle3, angle2, angle1] = cspice_m2eul(r, axis3, axis2, axis1)

   switch nargin
      case 4

         r     = zzmice_dp(r);
         axis3 = zzmice_int(axis3);
         axis2 = zzmice_int(axis2);
         axis1 = zzmice_int(axis1);

      otherwise

         error ( ['Usage: [_angle3_, _angle2_, _angle1_]  = '  ...
                  'cspice_m2eul( _r(3,3)_, axis3, axis2, axis1)']  )

   end

   %
   % Call the MEX library.
   %
   try
      [angle3, angle2, angle1] = mice('m2eul_c', r, axis3, axis2, axis1);
   catch
      rethrow(lasterror)
   end


