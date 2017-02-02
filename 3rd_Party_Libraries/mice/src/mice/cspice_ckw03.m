%-Abstract
%
%   CSPICE_CKW03 adds a type 3 segment to a C-kernel.
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
%      handle   file handle for an open CK file, returned from cspice_ckopn.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      begtim   encoded SCLK segment begin time.
%
%               [1,1] = size(begtim); double = class(begtim)
%
%      endtim   encoded SCLK segment end time.
%
%               [1,1] = size(endtim); double = class(endtim)
%
%      inst     NAIF instrument ID code.
%
%               [1,1] = size(inst); int32 = class(inst)
%
%      ref      name of the reference frame for the segment.
%
%               [1,c1] = size(ref), char = class(ref)
%
%                  or
%
%               [1,1] = size(ref), cell = class(ref)
%
%      avflag   a boolean signifying if the segment will contain
%               angular velocity.
%
%               [1,1] = size(avflag); logical = class(avflag)
%
%      segid    name to identify the segment.
%
%               [1,c2] = size(segid), char = class(segid)
%
%                  or
%
%               [1,1] = size(segid), cell = class(segid)
%
%      sclkdp   array containing the encoded SCLK times for the data.
%
%               [n,1] = size(endtim); double = class(endtim)
%
%      quats    array of SPICE style quaternions representing instrument
%               pointing.
%
%               [4,n] = size(endtim); double = class(endtim)
%
%      avvs     array of angular velocity vectors in units of radians per
%               second.
%
%               [3,n] = size(endtim); double = class(endtim)
%
%      starts   array containing the encoded SCLK interval start times of
%               each interpolation interval, the times must be strictly
%               increasing and coincide with pointing data times.
%
%               [m,1] = size(endtim); double = class(endtim)
%
%   the call:
%
%      cspice_ckw03( handle , ...
%                    begtime, ...
%                    endtime, ...
%                    inst   , ...
%                    ref    , ...
%                    avflag , ...
%                    segid  , ...
%                    sclkdp , ...
%                    quats  , ...
%                    avvs   , ...
%                    starts)
%
%   returns:
%
%      Adds a type 3 segment to a CK.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      INST3      = -77703;
%      NCOMCH     = 10;
%      REF        = 'J2000';
%      SEGID3     = 'Test type 3 test CK';
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
%      starts = [1:(MAXREC/2)]';
%      starts = (starts-1)*2*SPACING;
%
%      %
%      % Open a new CK, write the data, catch any errors.
%      %
%      try
%         handle = cspice_ckopn( 'test3.ck', 'ck', 0)
%         cspice_ckw03( handle , ...
%                       begtime, ...
%                       endtime, ...
%                       INST3  , ...
%                       REF    , ...
%                       avflag , ...
%                       SEGID3 , ...
%                       sclkdp , ...
%                       q      , ...
%                       expavvs, ...
%                       starts )
%      catch
%
%         error( [ 'Failure: ' lasterr] )
%      end
%
%      cspice_ckcls(handle)
%
%   MATLAB outputs:
%
%      The example code creates a CK with one type 3 segment.
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine ckw03_c.
%
%   MICE.REQ
%   CK.REQ
%   DAF.REQ
%   SCLK.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 11-JUL-2012, EDW (JPL)
%
%      Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 19-MAY-2006, EDW (JPL)
%
%-Index_Entries
%
%   write ck type_3 pointing data segment
%
%-&

function cspice_ckw03( handle , ...
                       begtime, ...
                       endtime, ...
                       inst   , ...
                       ref    , ...
                       avflag , ...
                       segid  , ...
                       sclkdp , ...
                       quats  , ...
                       avvs   , ...
                       starts )

   switch nargin
      case 11

         handle  = zzmice_int(handle);
         begtime = zzmice_dp(begtime);
         endtime = zzmice_dp(endtime);
         avflag  = zzmice_int(avflag);
         inst    = zzmice_int(inst);
         sclkdp  = zzmice_dp(sclkdp);
         quats   = zzmice_dp(quats);
         avvs    = zzmice_dp(avvs);
         starts  = zzmice_dp(starts);

      otherwise

         error ( [ 'Usage: '                ...
                   'cspice_ckw03( handle, ' ...
                           'begtime, '      ...
                           'endtime, '      ...
                           'inst, '         ...
                           'ref, '          ...
                           'avflag, '       ...
                           'segid, '        ...
                           'sclkdp(N), '    ...
                           'quats(4,N), '   ...
                           'avvs(3,N), '    ...
                           'starts(M))' ] )

   end

   %
   % Call the MEX library.
   %
   try
      mice( 'ckw03_c', handle, begtime,  endtime,  inst,  ref,  avflag,  ...
                       segid,  sclkdp,   quats,    avvs,  starts )
   catch
      rethrow(lasterror)
   end



