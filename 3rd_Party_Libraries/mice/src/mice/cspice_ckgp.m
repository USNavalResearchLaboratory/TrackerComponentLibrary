%-Abstract
%
%   CSPICE_CKGP returns pointing (attitude) for a specified
%   object at a user specified spacecraft clock time.
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
%                  The C-matrix returned by cspice_ckgp, if any, is the one
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
%      [ cmat, clkout, found ] = cspice_ckgp( inst, sclkdp, tol, ref )
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
%      clkout   the double precision scalar or double precision 1xN array of
%               encoded spacecraft clock time(s) associated with the returned
%               C-matrix 'cmat' (this value may differ from the requested
%               time, but never by more than the input tolerance 'tol')
%
%      found    the flag(s) indicating if the requested pointing is found.
%
%               [1,n] = size(found); logical = class(found)
%
%               'cmat', 'clkout' and 'found' return with the same vectorization
%               measure (N) as 'sclkdp'.
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
%      CK     = '/kernels/voyager2/ck/vg2_jup_qmw_na.bc';
%      cspice_furnsh( CK )
%
%      %
%      % Need to load a Voyager 2 SCLK kernel to convert from
%      % clock string to ticks.  Although not required for
%      % the Voyager spacecraft clocks, most modern spacecraft
%      % clocks require a leapseconds kernel to be loaded in
%      % addition to an SCLK kernel.
%      %
%      SCLK   = '/kernels/voyager2/sclk/vg200004.tsc';
%      cspice_furnsh( SCLK )
%
%      %
%      % The code for the Voyager 2 spacecraft clock is -32
%      %
%      SC     =  -32;
%
%      %
%      % The code for the narrow angle camera on the Voyager 2
%      % spacecraft is -32001.
%      %
%      INST   =  -32001;
%
%      %
%      % Spacecraft clock times for successive Voyager images always
%      % differ by more than 0:0:400.  This is an acceptable
%      % tolerance, and must be converted to "ticks" (units of
%      % encoded SCLK) for input to cspice_ckgp .
%      %
%      TOLVGR = '0:0:400';
%
%      %
%      % The reference frame we want is FK4.
%      %
%      REF    = 'FK4';
%
%      %
%      % The narrow angle camera boresight has direction
%      % ( 0, 0, 1 ) in camera (instrument) frame.
%      %
%      BORE   = [ 0.; 0.; 1.];
%      NPICS  = 2;
%
%      %
%      % Two Voyager 2 clock strings of interest.
%      %
%      SCLKCH =  strvcat( '2/18381:54:768' , '4/08970:58:768' );
%
%      %
%      % Convert tolerance from VGR formatted character string
%      % SCLK to ticks, which are units of encoded SCLK.
%      %
%      toltik = cspice_sctiks( SC, TOLVGR );
%
%      %
%      % cspice_ckgp requires encoded spacecraft clock time.
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
%      [ cmat, clkout, found ] = cspice_ckgp( INST, sclkdp, toltik, REF );
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
%            txt = sprintf( 'At SCLK time       : %s', SCLKCH(n,:) );
%            disp( txt )
%
%            txt = sprintf( 'The boresight vector: %16.8f %16.8f %16.8f', ...
%                                                                 bore_ref );
%            disp( txt )
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
%      At SCLK time       : 2/18381:54:768
%      The boresight vector:      -0.08938892       0.94794277       0.30563724
%
%      At SCLK time       : 4/08970:58:768 pointing not found
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ckgp_c.
%
%   MICE.REQ
%   CK.REQ
%   SCLK.REQ
%
%-Version
%
%   -Mice Version 1.2.0, 11-JUN-2013, EDW (JPL), SCK (JPL)
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
%   get ck pointing
%
%-&

function [cmat, clkout, found] = cspice_ckgp(inst, sclkdp, tol, ref)

   switch nargin
      case 4

         inst   = zzmice_int(inst);
         sclkdp = zzmice_dp(sclkdp);
         tol    = zzmice_dp(tol);
         ref    = zzmice_str(ref);

      otherwise

         error ( [ 'Usage: [_cmat(3,3)_, _clkout_, _found_] = ' ...
                   'cspice_ckgp(inst, _sclkdp_, tol, `ref`)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [cmat, clkout, found] = mice('ckgp_c', inst, sclkdp, tol, ref);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end

