%-Abstract
%
%   CSPICE_CKW01 adds a type 1 segment to a C-kernel.
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
%      handle   scalar integer value of the file handle for
%               an open CK file returned from cspice_ckopn
%
%      begtim   double precision scalar encoded SCLK segment
%               begin time
%
%      endtim   double precision scalar encoded SCLK segment
%               end time
%
%      inst     the scalar integer NAIF instrument ID code
%
%      ref      scalar string identifying the reference frame for the
%               segment
%
%      avflag   a scalar boolean indicating if the segment will contain
%               angular velocity
%
%      segid    a scalar string to identify the segment
%
%      sclkdp   double precision Nx1 array containing the encoded
%               SCLK times for the data
%
%      quats    a double precision 4xN matrix of SPICE style quaternions
%               representing instrument pointing
%
%      avvs     a double precision 3xN  matrix of angular
%               velocity vectors in units of radians per second
%
%   the call:
%
%      cspice_ckw01( handle , ...
%                    begtime, ...
%                    endtime, ...
%                    inst   , ...
%                    ref    , ...
%                    avflag , ...
%                    segid  , ...
%                    sclkdp , ...
%                    quats  , ...
%                    avvs )
%
%   returns:
%
%      Adds a type 1 segment to a CK.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      INST1      = -77701;
%      NCOMCH     = 10;
%      REF        = 'J2000';
%      SEGID1     = 'Test type 1 test CK';
%      SECPERTICK = 0.001;
%      SPACING    = 10.0;
%      MAXREC     = 50;
%
%      %
%      % Note, sclkdp is a vector input, not a vectorized scalar.
%      %
%      sclkdp    = [1:MAXREC]';
%      sclkdp    = (sclkdp - 1)*SPACING;
%
%      spinrate  = [1:MAXREC]*1.e-6;
%
%      theta     = [0:MAXREC-1]*SPACING;
%      theta     = theta .* spinrate;
%
%      %
%      % Create a zero-filled array for the angular velocity
%      % vectors. This allocates the needed memory and
%      % defines a variable of the correct shape.
%      %
%      expavvs = zeros( [3 MAXREC] );
%
%      a1 = zeros( [1 MAXREC] );
%      a2 = a1;
%
%      r  = cspice_eul2m( theta, a2, a1, 3, 1 ,3 );
%      q  = cspice_m2q( r );
%
%      %
%      % Fill the z component of the expavvs vectors with the
%      % corresponding spinrate element scaled to SECPERTICK.
%      %
%      expavvs(3,:) = spinrate/SECPERTICK;
%
%      begtime = sclkdp(1);
%      endtime = sclkdp(MAXREC);
%      avflag = 1;
%
%      %
%      % Open a new CK, write the data, catch any errors.
%      %
%      try
%         handle = cspice_ckopn( 'test1.ck', 'ck', 0)
%         cspice_ckw01( handle , ...
%                       begtime, ...
%                       endtime, ...
%                       INST1  , ...
%                       REF    , ...
%                       avflag , ...
%                       SEGID1 , ...
%                       sclkdp , ...
%                       q      , ...
%                       expavvs )
%      catch
%
%         error( [ 'Failure: ' lasterr] )
%      end
%
%      cspice_ckcls(handle)
%
%   MATLAB outputs:
%
%      The example code creates a CK with one type 1 segment.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ckw01_c.
%
%   MICE.REQ
%   CK.REQ
%   DAF.REQ
%   SCLK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   write ck type_1 pointing data segment
%
%-&

function cspice_ckw01( handle , ...
                       begtime, ...
                       endtime, ...
                       inst   , ...
                       ref    , ...
                       avflag , ...
                       segid  , ...
                       sclkdp , ...
                       quats  , ...
                       avvs )

   switch nargin
      case 10

         handle  = zzmice_int(handle);
         begtime = zzmice_dp(begtime);
         endtime = zzmice_dp(endtime);
         avflag  = zzmice_int(avflag);
         inst    = zzmice_int(inst);
         sclkdp  = zzmice_dp(sclkdp);
         quats   = zzmice_dp(quats);
         avvs    = zzmice_dp(avvs);

      otherwise

         error ( [ 'Usage: '                   ...
                   'cspice_ckw01( handle, '    ...
                                'begtime, '    ...
                                'endtime, '    ...
                                'inst, '       ...
                                'ref, '        ...
                                'avflag, '     ...
                                'segid, '      ...
                                'sclkdp(N), '  ...
                                'quats(4,N), ' ...
                                'avvs(3,N))' ] )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'ckw01_c', handle, begtime,  endtime,  inst,  ref,  ...
                       avflag, segid,    sclkdp,  quats,  avvs )
   catch
      rethrow(lasterror)
   end


