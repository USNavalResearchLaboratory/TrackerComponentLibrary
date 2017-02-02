%-Abstract
%
%   CSPICE_DELTET returns value of Delta ET (ET-UTC)
%   for an input epoch.
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
%      epoch    the double precision scalar or N-vector of
%               epochs at which "delta ET" is to be computed.
%               'epoch' may be either UTC or ephemeris seconds
%               past J2000, as specified by 'eptype'.
%
%      eptype   the scalar string indicating the type of input
%               epoch. It may be either of the following:
%
%                  'UTC'   UTC seconds past J2000 UTC.
%
%                  'ET'    Ephemeris seconds past J2000 TDB,
%                          also known as barycentric
%                          dynamical time (TDB).
%
%   the call:
%
%      delta = cspice_deltet( epoch, eptype )
%
%   returns:
%
%      delta   the double precision scalar or N-vector
%              values of
%
%                 "delta ET" = ET - UTC
%
%              at the input 'epoch'. This is added to UTC to
%              give ET, or subtracted from ET to give UTC.
%
%              'delta' return with the same vectorization measure
%              (N) as 'epoch'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load a leapsecond file.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % Define times of interest and the array size
%      % parameter.
%      %
%      SIZE     = 2004 - 1997 +1;
%      UTC_1997 = 'Jan 1 1997';
%      UTC_2004 = 'Jan 1 2004';
%
%      %
%      % Convert the UTC time strings to ET.
%      %
%      et_1997 = cspice_str2et( UTC_1997 );
%      et_2004 = cspice_str2et( UTC_2004 );
%
%      %
%      % Calculate the ET-UTC delta at Jan 1 1997
%      % and Jan 1 2004.
%      %
%      delt_1997 = cspice_deltet( et_1997, 'ET' );
%      delt_2004 = cspice_deltet( et_2004, 'ET' );
%
%      disp( 'Scalar:' )
%      disp( sprintf( 'Delta 1997: %f'  , delt_1997 ) )
%      disp( sprintf( 'Delta 2004: %f\n', delt_2004 ) )
%
%      %
%      % Given an array of 'SIZE' ephemeris times
%      % starting from value 'et_1997' with steps being
%      % of the number of seconds per Julian year, return
%      % the ET-UTC delta value for each time.
%      %
%      et   = [0:SIZE-1]*cspice_jyear + et_1997;
%      delt = cspice_deltet( et, 'ET' );
%
%      %
%      % Convert 'et' to 'utc'.
%      %
%      utc = cspice_et2utc( et, 'C', 3 );
%
%      disp( 'Vector:' )
%      for n=1:SIZE
%         txt = sprintf( '%s delta %f', utc(n,:), delt(n) );
%         disp( txt )
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Scalar:
%      Delta 1997: 62.183935
%      Delta 2004: 64.183912
%
%      Vector:
%      1997 JAN 01 00:00:00.000 delta 62.183935
%      1998 JAN 01 05:59:59.000 delta 63.183935
%      1999 JAN 01 11:59:58.000 delta 64.183935
%      2000 JAN 01 17:59:58.000 delta 64.183935
%      2000 DEC 31 23:59:58.000 delta 64.183934
%      2002 JAN 01 05:59:58.000 delta 64.183934
%      2003 JAN 01 11:59:58.000 delta 64.183934
%      2004 JAN 01 17:59:58.000 delta 64.183933
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine deltet_c.
%
%   MICE.REQ
%   TIME.REQ
%   KERNEL.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   difference between ephemeris time and utc
%
%-&

function [delta] = cspice_deltet( epoch, eptype )

   switch nargin
      case 2

         eptype = zzmice_str(eptype);
         epoch  = zzmice_dp(epoch);

      otherwise

         error ( 'Usage: [_delta_] = cspice_deltet( _epoch_, `eptype`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [delta] = mice('deltet_c', epoch, eptype) ;
   catch
      rethrow(lasterror)
   end


