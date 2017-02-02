%-Abstract
%
%   CSPICE_XF2EUL converts a state transformation matrix to
%   Euler angles and their derivatives with respect to a specified
%   set of axes.
%
%   The companion routine cspice_eul2xf converts Euler angles
%   and their derivatives with respect to a specified set of axes
%   to a state transformation matrix.
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
%      xform   a double precision 6x6 or double precision 6x6xN
%              array of a state transformation matrices from some frame
%              frame1 to another frame frame2
%
%      axisa
%      axisb
%      axisc   the scalar integers defining the axes desired for the
%              factorization of "r". All must be in the range from 1 to 3.
%              Moreover it must be the case that 'axisa' and 'axisb' are
%              distinct and that 'axisb' and 'axisc' are distinct.
%
%              Every rotation matrix can be represented as a product
%              of three rotation matrices about the principal axes
%              of a reference frame.
%
%                   r =  [ alpha ]     [ beta ]     [ gamma ]
%                                 axisa        axisb         axisc
%
%              The value 1 corresponds to the X axis.
%              The value 2 corresponds to the Y axis.
%              The value 3 corresponds to the Z axis.
%
%   the call:
%
%      [eulang, unique] = cspice_xf2eul(xform,  axisa, axisb, axisc)
%
%   returns:
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
%                The range of alpha and gamma is (-pi, pi].
%
%                The range of beta depends on the exact set of
%                axes used for the factorization.  For
%                factorizations in which the first and third axes
%                are the same, the range of beta is [0, pi].
%
%                For factorizations in which the first and third
%                axes are different, the range of beta is
%                [-pi/2, pi/2].
%
%                For rotations such that alpha and gamma are not
%                uniquely determined, alpha and dalpha/dt will
%                always be set to zero; gamma and dgamma/dt are
%                then uniquely determined.
%
%       unique   a boolean scalar or boolean 1XN array whether or not the
%                values in 'eulang' are uniquely determined.  If
%                the values are unique then 'unique' will be set to
%                true.  If the values are not unique and some
%                components ( eulang[1] and eulang[4] ) have value
%                zero, then 'unique' will have the value false.
%
%                'eulang' and 'unique' return with the same vectorization
%                measure (N) as 'xform'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%         KPL/MK
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            de421.bsp                     Planetary ephemeris
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( '/kernels/gen/lsk/naif0009.tls'
%                                '/kernels/gen/spk/de421.bsp'
%                                '/kernels/gen/pck/pck00009.tpc'
%                      )
%
%         \begintext
%
%
%      %
%      % Load the SPK, PCK and LSK kernels.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Suppose that you wish to determine the rate of change of
%      % the right ascension and declination of the pole of the moon,
%      % from the state transformation matrix that transforms J2000
%      % states to object fixed states.
%      %
%      % Using this routine with the routine sxform_c you can determine
%      % these instantaneous rates.
%      %
%      % Recall that the rotation component of tsipm is given by
%      %
%      %  [w]  [halfpi_c-dec] [ra+halfpi_c]
%      %     3               1             3
%      %
%      % Define the number of ephemeris times to perform the calculation.
%      %
%      N = 100;
%
%      %
%      % Calculate the separation of each ephemeris time, in seconds,
%      % over an eighteen year span.
%      %
%      STEP = 18 * 365 * cspice_spd/N;
%
%      %
%      % Base the ephemeris time set at May 15, 2007.
%      %
%      et = [0:N]*STEP +  cspice_str2et( 'May 15, 2007' );
%
%      %
%      % Calculate the state transformation matrices corresponding
%      % to 'et', then convert those matrices to Euler angles (3-1-3).
%      %
%      tsipm              = cspice_sxform( 'J2000', 'IAU_MOON', et );
%      [ eulang, unique ] = cspice_xf2eul( tsipm , 3, 1, 3 );
%
%      %
%      % From the Euler angles, calculate right ascension and declination.
%      % Form the UTC time string from 'et' (for output purposes).
%      %
%      ra  = eulang(3,:) - cspice_halfpi;
%      dec = cspice_halfpi - eulang(2,:);
%      utc = cspice_et2utc( et, 'c', 3 );
%
%      %
%      % As a convenience, output in a loop.
%      %
%      for m=1:N+1
%
%         if( unique(m) )
%            fprintf( 'UTC: %s\n', utc(m,:)                 )
%            fprintf( 'w        = %12.6f\n'  , eulang(1,m)  )
%            fprintf( 'dec      = %12.6f\n'  , dec(m)       )
%            fprintf( 'ra       = %12.6f\n'  , ra(m)        )
%            fprintf( 'd w/dt   = %14.9f\n'  , eulang(4,m)  )
%            fprintf( 'd dec/dt = %14.9f\n'  ,-eulang(5,m)  )
%            fprintf( 'd ra/dt  = %14.9f\n\n', eulang(6,m)  )
%         else
%            disp( 'The values in ''eulang'' not uniquely determined.' )
%         end
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%      %
%      % It is left as an exercise to the user to graphically display
%      % a plot of ra vs dec.
%      %
%      % plot(ra,dec)
%      %
%
%   MATLAB outputs:
%
%      The final five output sets, radians.
%
%               ...
%
%      UTC: 2024 AUG 20 04:48:00.002
%      w        =     2.719238
%      dec      =     1.188057
%      ra       =    -1.581646
%      d w/dt   =    0.000002658
%      d dec/dt =   -0.000000001
%      d ra/dt  =    0.000000004
%
%      UTC: 2024 OCT 24 21:36:00.003
%      w        =    -1.026611
%      dec      =     1.188773
%      ra       =    -1.576003
%      d w/dt   =    0.000002662
%      d dec/dt =    0.000000001
%      d ra/dt  =   -0.000000000
%
%      UTC: 2024 DEC 29 14:24:00.001
%      w        =     1.514167
%      dec      =     1.188984
%      ra       =    -1.573372
%      d w/dt   =    0.000002663
%      d dec/dt =   -0.000000001
%      d ra/dt  =   -0.000000001
%
%      UTC: 2025 MAR 05 07:12:00.000
%      w        =    -2.231256
%      dec      =     1.188290
%      ra       =    -1.567771
%      d w/dt   =    0.000002658
%      d dec/dt =    0.000000000
%      d ra/dt  =    0.000000004
%
%      UTC: 2025 MAY 10 00:00:00.000
%      w        =     0.307001
%      dec      =     1.188934
%      ra       =    -1.562882
%      d w/dt   =    0.000002663
%      d dec/dt =    0.000000001
%      d ra/dt  =   -0.000000001
%
%-Particulars
%
%   This function is intended to provide an inverse for the function
%   cspice_eul2xf.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine xf2eul_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 10-MAY-2011, EDW (JPL)
%
%      "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.1, 06-MAY-2009, EDW (JPL)
%
%      Added MICE.REQ reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 02-APR-2007, EDW (JPL)
%
%-Index_Entries
%
%   Euler angles and derivatives from state transformation
%
%-&

function [eulang, unique] = cspice_xf2eul(xform,  axisa, axisb, axisc)

   switch nargin
      case 4

         xform = zzmice_dp( xform);
         axisa = zzmice_int( axisa);
         axisb = zzmice_int( axisb);
         axisc = zzmice_int( axisc);

      otherwise

         error ( [ 'Usage: [_eulang(6)_, _unique_] = ' ...
                   'cspice_xf2eul(_xform(3,3)_, axisa, axisb, axisc)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [eulang, unique] = mice('xf2eul_c', xform, axisa, axisb, axisc);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      unique = zzmice_logical(unique);
   catch
      rethrow(lasterror)
   end


