%-Abstract
%
%   CSPICE_PXFORM returns the matrix that transforms position
%   vectors from one specified frame to another at a specified epoch.
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
%      from   the name of a reference frame in which a position is known.
%
%             [1,m] = size(from); char = class(from)
%
%      to     the name of a reference frame in which it is desired to represent
%             the position.
%
%             [1,l] = size(to); char = class(to)
%
%      et     epoch in ephemeris seconds past the epoch of J2000 (TDB) at which
%             the position transformation matrix should be evaluated.
%
%             [1,n] = size(et); double = class(et)
%
%   the call:
%
%      rotate = cspice_pxform( from, to, et )
%
%   returns:
%
%      rotate   operator(s) that transform position vector(s) from the
%               reference frame 'from' to frame 'to' at epoch 'et'
%
%               If [1,1] = size(et) then [3,3]   = size(rotate)
%               If [1,n] = size(et) then [3,3,n] = size(rotate)
%                                         double = class(rotate)
%
%               'rotate' returns with the same vectorization measure (N)
%               as 'et'.
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
%         File name: standard.tm
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
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0009.tls'
%                                'de421.bsp'
%                                'pck00009.tpc' )
%
%         \begintext
%
%      %
%      % Output the right ascension and declination of the earth's pole
%      % in the J2000 frame approximately every month for the time
%      % interval January 1, 1990 to January 1, 2010 (UTC).
%      %
%      %
%      % Load a standard kernel set.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define the time bounds for the time interval,
%      % 20 years, convert to ephemeris time J2000.
%      %
%      utc_bounds = strvcat( '1 Jan 1990', '1 Jan 2010' );
%      et_bounds  = cspice_str2et( utc_bounds );
%
%      %
%      % Step in units of a month. 20 years ~ 240 months.
%      %
%      step = (et_bounds(2) - et_bounds(1) ) / 240.;
%
%      %
%      % Create an array of 240 ephemeris times ending at
%      % ~et_bound(2) in intervals of 'step'.
%      %
%      et = [1:240]*step + et_bounds(1);
%
%      %
%      % Set the conversion constant "radians to degrees."
%      %
%      r2d = cspice_dpr;
%
%      %
%      % Convert the 240-vector of 'et' to an array of corresponding
%      % transformation matrices (dimensions (3,3,240) ).
%      %
%      mat = cspice_pxform( 'IAU_EARTH', 'J2000', et );
%
%      %
%      % Extract the pole vector from the transformation matrix,
%      % convert to RA and DEC expressed in degrees.
%      %
%
%      %
%      % The last column in each matrix is the pole vector (z = (0,0,1))
%      % of the earth in IAU expressed in J2000.
%      %
%      % Recall, MATLAB uses 1 based indexing, so (:,3,:) represents.
%      % the third column of the matrices.
%      %
%      pole = mat(:,3,:);
%
%      %
%      % 'pole' ready for use in cspice_radrec.
%      %
%      [radius, ra, dec] = cspice_recrad( pole );
%
%      %
%      % Output the ephemeris time and the corresponding
%      % angular values (in degrees). 'ra' and 'dec' return
%      % as double precision 240-vectors.
%      %
%      ra  = ra  * r2d;
%      dec = dec * r2d;
%
%      %
%      % Create an array of values for output.
%      %
%      output = [  et; ra; dec ];
%
%      fprintf( '%24.8f %16.8f %16.8f\n', output );
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      A partial output centered on et = 0:
%
%                              ...
%
%      -18408539.52023917     180.00373915      89.99675084
%      -15778739.49107254     180.00320499      89.99721501
%      -13148939.46190590     180.00267082      89.99767918
%      -10519139.43273926     180.00213665      89.99814334
%       -7889339.40357262     180.00160249      89.99860751
%       -5259539.37440598     180.00106832      89.99907168
%       -2629739.34523934     180.00053415      89.99953584
%             60.68392730     359.99999999      89.99999999
%        2629860.71309394     359.99946582      89.99953582
%        5259660.74226063     359.99893165      89.99907166
%        7889460.77142727     359.99839749      89.99860749
%       10519260.80059391     359.99786332      89.99814332
%       13149060.82976055     359.99732915      89.99767916
%       15778860.85892719     359.99679499      89.99721499
%       18408660.88809383     359.99626082      89.99675082
%
%                              ...
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pxform_c.
%
%   MICE.REQ
%   ROTATION.REQ
%   FRAMES.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 09-NOV-2012, EDW (JPL), SCK (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   Find a position transformation matrix
%
%-&

function [rotate] = cspice_pxform(from, to, et)

   switch nargin
      case 3

         from = zzmice_str(from);
         to   = zzmice_str(to);
         et   = zzmice_dp(et);

      otherwise

         error ( [ 'Usage: [_rotate(3,3)_] = ' ...
                   'cspice_pxform( `from`, `to`, _et_ )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [rotate] = mice('pxform_c',from,to,et);
   catch
      rethrow(lasterror)
   end




