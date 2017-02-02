%-Abstract
%
%   CSPICE_CKGPAV returns pointing (attitude) and angular velocity
%   for a specified object at a user specified spacecraft clock time.
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
%      inst     NAIF ID for the instrument, spacecraft, or other structure for
%               which pointing is requested. The frame fixed to this object is
%               called the "instrument frame" or "instrument-fixed" frame.
%
%               [1,1] = size(inst); int32 = class(inst)
%
%      sclkdp   encoded spacecraft clock time(s) for which pointing is
%               requested.
%
%               [1,n] = size(sclkdp); double = class(sclkdp)
%
%      tol      time tolerance given in ticks (+/-), the units of encoded
%               spacecraft clock time, about 'sclkdp'
%
%               [1,1] = size(tol); double = class(tol)
%
%                  The C-matrix returned by cspice_ckgpav, if any, is the one
%                  whose time tag is closest to 'sclkdp' and within 'tol'
%                  units of 'sclkdp'.
%
%               In general, because using a non-zero tolerance affects
%               selection of the segment from which the data is obtained,
%               users are strongly discouraged from using a non-zero
%               tolerance when reading CKs with continuous data. Using
%               a non-zero tolerance should be reserved exclusively to
%               reading CKs with discrete data because in practice
%               obtaining data from such CKs using a zero tolerance is
%               often not possible due to time round off.
%
%      ref      naming the desired reference frame for the returned pointing.
%
%               [1,m] = size(ref); char = class(ref)
%
%   the call:
%
%      [ cmat, av, clkout, found ] = cspice_ckgpav( inst, sclkdp, tol, ref)
%
%   returns:
%
%      cmat     rotation matrix(ces) that transform components of a vector
%               expressed in the frame specified by 'ref' to components
%               expressed in the frame tied to the instrument, spacecraft, or
%               other structure at time(s) 'clkout'.
%
%               If [1,1] = size(sclkdp) then [3,3]   = size(cmat)
%               If [1,n] = size(sclkdp) then [3,3,n] = size(cmat)
%                                             double = class(cmat)
%
%      av       angular velocity measured in radians per second (this is the
%               axis about which the reference frame tied to the instrument is
%               rotating in the right-handed sense at time 'clkout')
%
%               [3,n] = size(av); double = class(av)
%
%      clkout   encoded spacecraft clock time(s) associated with the returned
%               C-matrix 'cmat' (this value may differ from the requested time,
%               but never by more than the input tolerance 'tol').
%
%               [1,n] = size(clkout); double = class(clkout)
%
%      found    the flag(s) indicating if the requested pointing is found.
%
%               [1,n] = size(found); logical = class(found)
%
%               'cmat', 'av', 'clkout', and 'found' return with the same
%               vectorization measure (N) as 'sclkdp'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load the CK pointing file.
%      %
%      CK     = '/kernels/cassini/ck/030810_031019_c39_port1_pa.bc';
%      cspice_furnsh( CK )
%
%      %
%      % Need to load a Cassini SCLK kernel to convert from
%      % clock string to ticks and the leapseconds kernel.
%      %
%      SCLK   = '/kernels/cassini/sclk/cas00082.tsc';
%      cspice_furnsh( SCLK )
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      % The code for the Cassini spacecraft is -82
%      %
%      SC     = -82;
%
%      % The code for the Cassini spacecraft bus is -82000.
%      %
%      INST   = -82000;
%
%      %
%      % The reference frame we want is J2000.
%      %
%      REF    = 'J2000';
%
%      %
%      % The bus boresight has direction ( 0, 0, 1 ) in the
%      % bus (instrument) frame.
%      %
%      BORE   = [ 0.; 0.; 1.];
%
%      %
%      % The spacecraft clock ticks every 1/256 of a second.
%      %
%      TICS   =  256.;
%
%      %
%      % Two Cassini clock strings of interest.
%      %
%      SCLKCH =  strvcat( '1440000000.000' , '1445000000.000' );
%
%      %
%      % Calculate a tolerance for pointing data. Use 5 seconds,
%      % convert to ticks by multiplying by the number of ticks
%      % per second.
%      %
%      toltik = 5. * TICS;
%
%      %
%      % cspice_ckgpav requires encoded spacecraft clock time.
%      %
%      sclkdp = cspice_scencd( SC, SCLKCH );
%
%      %
%      % Retrieve the 'REF' reference frame to 'INST' reference frame
%      % transformation matrix at time sclkdp with a tolerance
%      % 'toltik'.
%      %
%      %   [INST] = [cmat][ref]
%      %
%      [ cmat, av, clkout, found ] = cspice_ckgpav( INST, sclkdp, toltik, REF);
%
%      for n=1:2
%
%         if( found(n) )
%
%            %
%            % Transform the 'BORE' vector from 'INST' reference frame to
%            % 'REF' frame.
%            %                T
%            %  [ref] = [cmat] [INST]
%            %
%            bore_ref = cmat(:,:,n)' * BORE;
%
%            txt = sprintf( 'At SCLK time               : %s', SCLKCH(n,:) );
%            disp( txt )
%
%            txt = sprintf(                                             ...
%                  'The boresight vector       : %11.8f %11.8f %11.8f', ...
%                                                               bore_ref );
%            disp( txt )
%
%            txt = sprintf(                                            ...
%                 'The angular velocity vector: %10.7e %10.7e %10.7e', ...
%                                                               av(:,n) );
%            disp( txt )
%
%            disp( ' ' )
%
%         else
%
%            txt = sprintf( 'At SCLK time       : %s pointing not found', ...
%                                                              SCLKCH(n,:) );
%            disp( txt )
%            disp( ' ' )
%
%         end
%
%      end
%
%   MATLAB outputs:
%
%      At SCLK time               : 1440000000.000
%      The boresight vector       : -0.25517775  0.88933091  0.37944018
%      The angular velocity vector: -3.1532789e-13 -7.8841788e-09 1.7857623e-08
%
%      At SCLK time               : 1445000000.000
%      The boresight vector       : -0.31149356  0.87452887  0.37171362
%      The angular velocity vector: -8.0798136e-11 -7.1102943e-10 1.5325493e-09
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ckgpav_c.
%
%   MICE.REQ
%   CK.REQ
%   SCLK.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 08-NOV-2012, EDW (JPL), SCK (JPL)
%
%       I/O descriptions edits to conform to Mice documentation format.
%
%      "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.1.1, 03-JUN-2010, BVS (JPL)
%
%      Edits to header. Added warning regarding non-zero tolerance.
%
%   -Mice Version 1.1.0, 23-FEB-2009, EDW (JPL)
%
%      Added zzmice_str call on input 'ref' to convert string cells to
%      character arrays if 'ref' has type string cells. Added proper
%      markers for usage string variable types.
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   get ck pointing and angular velocity
%
%-&

function [cmat, av, clkout, found] = cspice_ckgpav(inst, sclkdp, tol, ref)

   switch nargin
      case 4

         inst   = zzmice_int(inst);
         sclkdp = zzmice_dp(sclkdp);
         tol    = zzmice_dp(tol);
         ref    = zzmice_str(ref);

      otherwise

         error ( [ 'Usage: [_cmat(3,3)_, _av(3)_, _clkout_, _found_] = ' ...
                   'cspice_ckgpav(inst, _sclkdp_, tol, `ref`)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [cmat, av, clkout, found] = mice('ckgpav_c', inst, sclkdp, tol, ref);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end



