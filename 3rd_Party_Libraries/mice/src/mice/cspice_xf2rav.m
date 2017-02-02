%-Abstract
%
%   CSPICE_XF2RAV determines the rotation matrix and angular
%   velocity of the rotation from a state transformation matrix.
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
%              array of a state transformations from some frame
%              "frame1" to another frame "frame2"
%
%   the call:
%
%      [rot, av] = cspice_xf2rav(xform)
%
%   returns:
%
%      rot   a double precision 3x3 or double precision 3x3xN array of
%            rotation matrices that gives the transformation from some
%            frame "frame1" to another frame "frame2"
%
%      av    the double precision 3x1 array or double precision
%            3xN array of angular velocities of the transformation
%
%            If 'p' is the position of a fixed point in "frame2,"
%            then from the point of view of "frame1," 'p' rotates
%            (in a right handed sense) about an axis parallel to
%            'av'.  Moreover the rate of rotation in radians per unit
%            time is given by the length of 'av'.
%
%            More formally, the velocity 'v' of 'p' in "frame1" is
%            given by
%                                  t
%               v  = av x ( rot * p )
%
%            The components of 'av' are given relative to "frame1."
%
%            'rot' and 'av' return with the same vectorization
%            measure (N) as 'xform'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      %  Load a set of kernels: an SPK file, a PCK file
%      %  and a leapseconds file. Use a meta kernel
%      %  for convenience.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Create an array of 10001 ephemeris times based at July 1 2007.
%      %
%      et    = [0: 10000]* cspice_spd + cspice_str2et( 'July 1 2007' );
%
%      %
%      % Calculate the state transformation matrices from J2000 to IAU_MOON
%      % for 'et'.
%      %
%      xform = cspice_sxform( 'J2000', 'IAU_MOON', et );
%
%      %
%      % Convert the set of 'xform' matrices to the corresponding rotation
%      % matrices and angular velocity vectors.
%      %
%      [ rot, av ] = cspice_xf2rav(xform);
%
%      %
%      % Use the converted outputs from cspice_xf2rav to recompute a set
%      % of state transformation matrices.
%      %
%      strans = cspice_rav2xf( rot, av );
%
%      %
%      % Calculate the maximum value of the absolute difference between
%      % 'xform' and 'strans'.
%      %
%      disp( 'Vector:' )
%      fprintf(                                                              ...
%         'Maximum absolute difference between rotation matrices: %8.6e\n', ...
%                                   max( max( max( abs(strans - xform) ) ) )   )
%
%      %
%      %  It's always good form to unload kernels after use,
%      %  particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Maximum absolute difference between rotation matrices: 1.694066e-21
%
%      In this case, a value on the order of -21 indicates numerical
%      equivalence.
%
%-Particulars
%
%   This routine is an inverse of the routine cspice_rav2xf.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine xf2rav_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 06-MAY-2009, EDW (JPL)
%
%      Added MICE.REQ reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 11-APR-2007, EDW (JPL)
%
%-Index_Entries
%
%   State transformation to rotation and angular velocity
%
%-&

function [rot, av] = cspice_xf2rav(xform)

   switch nargin
      case 1

         xform = zzmice_dp(xform);

      otherwise

         error ( ['Usage: [_rot(3,3)_, _av(3)_] = ' ...
                  'cspice_xf2rav(_xform(6,6)_)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [rot, av] = mice('xf2rav_c', xform);
   catch
      rethrow(lasterror)
   end


