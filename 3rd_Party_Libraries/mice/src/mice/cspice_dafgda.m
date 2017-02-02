%-Abstract
%
%   CSPICE_DAFGDA reads the double precision data bounded by two addresses
%   within a DAF.
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
%      handle   file handle referring to a DAF.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      baddr,
%      eaddr    initial and final addresses of a contiguous set of double
%               precision numbers within a DAF. Presumably, these make up
%               all or part of a particular array.
%
%               Note that DAF addresses begin at 1 as in the
%               FORTRAN version of the SPICE Toolkit.
%
%               [1,1] = size(baddr); int32 = class(baddr)
%               [1,1] = size(eaddr); int32 = class(eaddr)
%
%   the call:
%
%      data = cspice_dafgda( handle, baddr, eaddr )
%
%   returns:
%
%      data   are the double precision data contained between
%             the specified addresses within the specified file.
%
%             'data' has length = end - begin + 1.
%
%             [1,length] = size(data); double = class(data)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Open the type 8 SPK "gda.bsp" for read access then read the
%      % data from the first segment. The segment contains 100
%      % 6 element records plus four additional elements.
%      %
%      handle = cspice_dafopr( 'gda.bsp');
%
%      %
%      % Begin a forward search; find the first segment; read the
%      % segment summary.
%      %
%      cspice_dafbfs( handle )
%      found    = cspice_daffna;
%      [dc, ic] = cspice_dafgs( 2, 6 );
%
%      %
%      % Retrieve the data begin and end addresses.
%      %
%      baddr = ic(5);
%      eaddr = ic(6);
%
%      fprintf( 'Beginning address       : %d\n', baddr )
%      fprintf( 'Ending address          : %d\n', eaddr )
%      fprintf( 'Number of data elements : %d\n', eaddr - baddr + 1 )
%
%      %
%      % Extract all data bounded by the begin and end addresses.
%      %
%      data = cspice_dafgda( handle, baddr, eaddr );
%
%      %
%      % Check 'data'. It should show an array of 604 doubles (4 + 6 * 100).
%      %
%      fprintf( 'Size of data array      : ' )
%      fprintf( '%d ', size(data) )
%      fprintf('\n\n')
%
%      %
%      % Check the data. Each set of 6 element records should possess the
%      % property:
%      %
%      %   record(6) = record(6)  + 1000.
%      %        i            i-1
%      %
%      fprintf( ' %7.2f ', data(1:6) )
%      fprintf('\n')
%
%      fprintf( ' %6.2f ', data(7:12) )
%      fprintf('\n')
%
%      %
%      % SAFELY close the file
%      %
%      cspice_dafcls(handle)
%
%   Matlab outputs:
%
%      Beginning address       : 385
%      Ending address          : 988
%      Number of data elements : 604
%      Size of data array      : 1 604
%
%          0.00     1.00     2.00     3.00     4.00     5.00
%       1000.00  1001.00  1002.00  1003.00  1004.00  1005.00
%
%   cspice_dafgda returned 604 double precision data values between DAF
%   addresses 385 and 988. The second 6-vector shows the property of 1000
%   more than the previous set, as expected.
%
%-Particulars
%
%   The principal reason that DAFs are so easy to use is that
%   the data in each DAF are considered to be one long contiguous
%   set of double precision numbers. You can grab data from anywhere
%   within a DAF without knowing (or caring) about the physical
%   records in which they are stored.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine dafgda_c.
%
%   MICE.REQ
%   DAF.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 17-JUL-2012, EDW (JPL)
%
%-Index_Entries
%
%   read data from DAF address
%
%-&

function [data] = cspice_dafgda( handle, baddr, eaddr)

   switch nargin
      case 3

         handle = zzmice_int(handle);
         baddr  = zzmice_int(baddr);
         eaddr  = zzmice_int(eaddr);

      otherwise

         error ( 'Usage: data = cspice_dafgda( handle, baddr, eaddr)' )

   end

   %
   % Call the MEX library.
   %
   try
      [data] = mice( 'dafgda_c', handle, baddr, eaddr );
   catch
      rethrow(lasterror)
   end

