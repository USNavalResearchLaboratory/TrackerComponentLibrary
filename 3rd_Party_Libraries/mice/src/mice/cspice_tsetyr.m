%-Abstract
%
%   CSPICE_TSETYR sets the lower bound on the 100 year range.
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
%      year   the scalar integer year associated with the lower
%             bound on all year expansions
%
%   the call:
%
%      cspice_tsetyr( year )
%
%   returns:
%
%      None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Set the lower bound for two digit century representation to 1980.
%
%      cspice_tsetyr(1980)
%
%   This cause the SPICE time system to interpret two digit year
%   values in the following manner:
%
%      year value      implied value
%      ----------      -----------
%      00              2000
%      01              2001
%       .                .
%       .                .
%      79              2079  <- lower bound + 99 years
%      80              1980  <- lower bound
%       .                .
%       .                .
%       .                .
%      99              1999
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine tsetyr_c.
%
%   MICE.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 12-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   Set the interval of expansion for abbreviated years
%
%-&

function cspice_tsetyr( year )

  switch nargin
      case 1

         year = zzmice_int(year);

      otherwise

         error ( 'Usage: cspice_tsetyr(year)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('tsetyr_c', year);
   catch
      rethrow(lasterror)
   end


