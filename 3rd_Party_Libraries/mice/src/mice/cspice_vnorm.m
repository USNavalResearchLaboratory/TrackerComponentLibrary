
%-Abstract
%
%   CSPICE_VNORM returns the magnitude of a double precision, 3-dimensional
%   array or set of such arrays.
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
%      v1   any arbitrary 3-vector(s).
%
%           [3,n] = size(v1); double = class(v1)
%
%   the call:
%
%      vnorm = cspice_vnorm(v1)
%
%   returns:
%
%      vnorm   the positive definite magnitude(s) of 'v1', ||v1||, calculated
%              in a numerically stable way.
%
%              [1,n] = size(vnorm); double = class(vnorm)
%
%              'vnorm' returns with the same measure of vectorization (N)
%              as 'v1'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Set a size for the Nx1 array of ephemeris times.
%      %
%      N = 1000;
%
%      %
%      %  Load a set of kernels: an SPK file, a PCK
%      %  file and a leapseconds file. Use a meta
%      %  kernel for convenience.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Set a reference epoch, convert the string representation
%      % to ET.
%      %
%      utc = 'Jan 1 2010 12:34:56';
%      et_0 = cspice_str2et( utc );
%
%      %
%      % Create an array of N elements off the reference epoch in
%      % steps of one day in ET seconds.
%      %
%      et = [1:N]*cspice_spd() + et_0;
%
%      %
%      % Calculate the geometric position of Mercury with respect to
%      % the earth, without aberration correction, at time 'et'.
%      %
%      target   = 'Mercury';
%      frame    = 'J2000';
%      abcorr   = 'none';
%      observer = 'Earth';
%
%      disp('Scalar')
%
%      [pos, ltime] = cspice_spkpos( target, et_0, frame, abcorr, observer );
%
%      %
%      % Calculate the  magnitude of the position vector returned
%      % from cspice_spkpos.
%      %
%      vmag = [ cspice_vnorm( pos ) ]'
%
%
%      disp('Vectorized')
%
%      [pos, ltime] = cspice_spkpos( target, et, frame, abcorr, observer );
%
%      %
%      % Calculate the 1xN array of magnitudes of the N position vectors
%      % returned from cspice_spkpos.
%      %
%      vmag = [ cspice_vnorm( pos ) ]'
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Matlab due to data persistence.
%      %
%      cspice_kclear
%
%   Matlab outputs:
%
%      Scalar
%
%      vmag =
%
%           1.041179020016839e+08
%
%      Vectorized
%
%      vmag =
%
%           1.025277808532095e+08
%           1.013810734850024e+08
%           1.006917331078249e+08
%           1.004611905399685e+08
%           1.006785378087216e+08
%           1.013217510356714e+08
%
%                   ...
%
%           2.076962759877405e+08
%           2.072238879729207e+08
%           2.066729526239417e+08
%           2.060458187524104e+08
%           2.053445529777324e+08
%           2.045709570882183e+08
%
%-Particulars
%
%   The magnitude calculation takes care to avoid overflow while computing
%   the norm of the input vector 'v1'. The logic determines the component of
%   'v1' whose magnitude is the largest. Calling this magnitude v1max, the
%   norm is computed using the formula
%
%       vnorm  =  v1max *  ||  (1/v1max) * v1  ||
%
%   where the notation ||x|| indicates the norm of the vector x.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine vnorm_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 24-APR-2010 (EDW)
%
%-Index_Entries
%
%   norm of 3-dimensional vector
%
%-&

function [vnorm] = cspice_vnorm(v1)

   switch nargin
      case 1

         v1 = zzmice_dp(v1);

      otherwise

         error ( 'Usage: [_vnorm_] = cspice_vnorm(_v1(3)_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [vnorm] = mice('vnorm_c', v1);
   catch
      rethrow(lasterror)
   end
