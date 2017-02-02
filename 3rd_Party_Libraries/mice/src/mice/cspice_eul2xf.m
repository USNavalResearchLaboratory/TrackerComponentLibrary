%-Abstract
%
%   CSPICE_EUL2XF computes a state transformation from an Euler angle
%   factorization of a rotation and the derivatives of those Euler
%   angles.
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
%       eulang   the double precision 6x1 or double precision 6xN
%                array of Euler angles corresponding to the
%                specified factorization
%
%                If we represent r as shown here:
%
%                   r =  [ alpha ]     [ beta ]     [ gamma ]
%                                 axisa        axisb         axisc
%
%                then (6x1)
%
%                  eulang[1] = alpha
%                  eulang[2] = beta
%                  eulang[3] = gamma
%                  eulang[4] = dalpha/dt
%                  eulang[5] = dbeta/dt
%                  eulang[6] = dgamma/dt
%
%                or (6xN)
%
%                  eulang[:,N] = alpha_N
%                  eulang[:,N] = beta_N
%                  eulang[:,N] = gamma_N
%                  eulang[:,N] = dalpha_N/dt
%                  eulang[:,N] = dbeta_N/dt
%                  eulang[:,N] = dgamma_N/dt
%
%      axisa
%      axisb
%      axisc     the scalar integers defining the axes desired for the
%                factorization of "r". All must be in the range from 1 to 3.
%
%                Every rotation matrix can be represented as a product
%                of three rotation matrices about the principal axes
%                of a reference frame.
%
%                   r =  [ alpha ]     [ beta ]     [ gamma ]
%                                 axisa        axisb         axisc
%
%                The value 1 corresponds to the X axis.
%                The value 2 corresponds to the Y axis.
%                The value 3 corresponds to the Z axis.
%
%   the call:
%
%      xform = cspice_eul2xf(eulang, axisa, axisb, axisc)
%
%   returns:
%
%      xform   a double precision 6x6 or double precision 6x6xN
%              array of a state transformation matrices corresponding to
%              'r' as described above
%
%              'xform'  returns with the same vectorization
%              measure (N) as 'eulang'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Suppose you have a set of Euler angles and their derivatives
%      % for a 3 1 3 rotation, and that you would like to determine
%      % the equivalent angles and derivatives for a 1 2 3 rotation.
%      %
%      % r = [alpha]  [beta]  [gamma]
%      %           3       1        3
%      %
%      % r = [roll]  [pitch]  [yaw]
%      %           1        2      3
%      %
%      % The following code fragment will perform the desired computation.
%      %
%      abgang = [0.01; 0.03; 0.09; -0.001; -0.003; -0.009 ];
%
%      xform              = cspice_eul2xf( abgang, 3, 1, 3 );
%      [ rpyang, unique ] = cspice_xf2eul( xform , 1, 2, 3 );
%
%      if( unique )
%         disp( '1-2-3 equivalent rotation to input (radians):')
%         fprintf( 'Roll  %12.6f, dRoll/dt  %12.6f\n', rpyang(1), rpyang(4) )
%         fprintf( 'Pitch %12.6f, dPitch/dt %12.6f\n', rpyang(2), rpyang(5) )
%         fprintf( 'Yaw   %12.6f, dYaw/dt   %12.6f\n', rpyang(3), rpyang(6) )
%      else
%         disp( 'The values in ''rpyang'' not uniquely determined.' )
%      end
%
%   MATLAB outputs:
%
%      1-2-3 equivalent rotation to input (radians):
%      Roll      0.029999, dRoll/dt     -0.003000
%      Pitch    -0.000300, dPitch/dt     0.000060
%      Yaw       0.099996, dYaw/dt      -0.009999
%
%-Particulars
%
%   This function is intended to provide an inverse for the function
%   cspice_xf2eul.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine eul2xf_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.2, 29-FEB-2012, EDW (JPL)
%
%      Edit to "Usage" string. "xform(3,3)" corrected to read
%      "xform(6,6)."
%
%   -Mice Version 1.0.1, 06-MAY-2009, EDW (JPL)
%
%      Added MICE.REQ reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 02-APR-2007, EDW (JPL)
%
%-Index_Entries
%
%   State transformation from Euler angles and derivatives
%
%-&

function [xform] = cspice_eul2xf(eulang, axisa, axisb, axisc)

   switch nargin
      case 4

         eulang= zzmice_dp(eulang);
         axisa = zzmice_int( axisa);
         axisb = zzmice_int( axisb);
         axisc = zzmice_int( axisc);

      otherwise

         error ( ['Usage: [_xform(6,6)_] = ' ...
                 'cspice_eul2xf(_eulang(6)_, axisa, axisb, axisc)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [xform] = mice('eul2xf_c', eulang, axisa, axisb, axisc);
   catch
      rethrow(lasterror)
   end



