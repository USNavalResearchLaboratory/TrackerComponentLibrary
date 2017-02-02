%-Abstract
%
%   CSPICE_SPKPOS returns the position of a target body relative
%   to an observing body, optionally corrected for light time
%   (planetary aberration) and stellar aberration.
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
%      targ      the scalar string name of a target body.
%                Optionally, you may supply the integer ID code
%                for the object as an integer string, i.e. both
%                'MOON' and '301' are legitimate strings that
%                indicate the Moon is the target body.
%
%                The target and observer define a position vector
%                whose position component points from the observer
%                to the target.
%
%                [1,c1] = size(target), char = class(target)
%
%                   or
%
%                [1,1] = size(target); cell = class(target)
%
%      et        the scalar or 1XN-vector of double precision ephemeris
%                time, expressed as seconds past J2000 TDB, at which
%                position of the target body relative to the observer
%                the is to be  computed, 'et' refers to time at the
%                observer's location.
%
%                [1,N] = size(et), double = class(et)
%
%      ref       the scalar string name of the reference frame relative
%                to which the output position vector should be
%                expressed. This may be any frame supported by the SPICE
%                system, including built-in frames (documented in the
%                Frames Required Reading) and frames defined by a loaded
%                frame kernel (FK).
%
%                When 'ref' designates a non-inertial frame, the
%                orientation of the frame is evaluated at an epoch
%                dependent on the selected aberration correction.
%
%                [1,c2] = size(ref), char = class(ref)
%
%                   or
%
%                [1,1] = size(ref); cell = class(ref)
%
%      abcorr    a scalar string that indicates the aberration corrections
%                to apply to the position of the target body to account
%                for one-way light time and stellar aberration.
%
%                'abcorr' may be any of the following:
%
%                   'NONE'     Apply no correction. Return the
%                              geometric position of the target
%                              body relative to the observer.
%
%                The following values of 'abcorr' apply to the
%                "reception" case in which photons depart from the
%                target's location at the light-time corrected epoch
%                et-lt and *arrive* at the observer's location at
%                'et':
%
%                   'LT'       Correct for one-way light time (also
%                              called "planetary aberration") using a
%                              Newtonian formulation. This correction
%                              yields the position of the target at the
%                              moment it emitted photons arriving at
%                              the observer at 'et'.
%
%                              The light time correction uses an
%                              iterative solution of the light time
%                              equation (see Particulars for details).
%                              The solution invoked by the "LT" option
%                              uses one iteration.
%
%                   'LT+S'     Correct for one-way light time and
%                              stellar aberration using a Newtonian
%                              formulation. This option modifies the
%                              position obtained with the "LT" option to
%                              account for the observer's velocity
%                              relative to the solar system
%                              barycenter. The result is the apparent
%                              position of the target---the position
%                              of the target as seen by the
%                              observer.
%
%                   'CN'       Converged Newtonian light time correction.
%                              This option produces a solution that is at
%                              least as accurate at that obtainable 
%                              with the 'LT' option. Whether the 'CN' 
%                              solution is substantially more accurate 
%                              depends on the geometry of the 
%                              participating objects and on the 
%                              accuracy of the input data. In all 
%                              cases this routine will execute more 
%                              slowly when a converged solution is 
%                              computed. See the section titled "The 
%                              Computation of Light Time" in the SPK
%                              Required Reading document spk.req for 
%                              details.
%
%                   'CN+S'     Converged Newtonian light time
%                              correction and stellar aberration
%                              correction.
%
%
%                The following values of 'abcorr' apply to the
%                "transmission" case in which photons *depart* from
%                the observer's location at 'et' and arrive at the
%                target's location at the light-time corrected epoch
%                et+lt:
%
%                   'XLT'      "Transmission" case:  correct for
%                              one-way light time using a Newtonian
%                              formulation. This correction yields the
%                              position of the target at the moment it
%                              receives photons emitted from the
%                              observer's location at 'et'.
%
%                   'XLT+S'    "Transmission" case:  correct for
%                              one-way light time and stellar
%                              aberration using a Newtonian
%                              formulation  This option modifies the
%                              position obtained with the "XLT" option to
%                              account for the observer's velocity
%                              relative to the solar system
%                              barycenter. The position indicates the
%                              direction that photons emitted from the
%                              observer's location must be "aimed" to
%                              hit the target.
%
%                   'XCN'      "Transmission" case:  converged
%                              Newtonian light time correction.
%
%                   'XCN+S'    "Transmission" case:  converged
%                              Newtonian light time correction and
%                              stellar aberration correction.
%
%
%                Neither special nor general relativistic effects are
%                accounted for in the aberration corrections applied
%                by this routine.
%
%                Neither letter case or embedded blanks are significant
%                in the 'abcorr' string.
%
%                [1,c3] = size(abcorr), char = class(abcorr)
%
%                   or
%
%                [1,1] = size(abcorr); cell = class(abcorr)
%
%      obs       the scalar string name of a observing body.
%                Optionally, you may supply the integer ID code
%                for the object as an integer string, i.e. both
%                'MOON' and '301' are legitimate strings that
%                indicate the Moon is the observing body.
%
%                [1,c4] = size(obs), char = class(obs)
%
%                   or
%
%                [1,1] = size(obs); cell = class(obs)
%
%   the call:
%
%      [pos, lt] = cspice_spkpos(targ, et, ref, abcorr, obs)
%
%   returns:
%
%      pos   a double precision 3x1 array or double precision
%            3xN array representing the position of the target body in
%            kilometers relative to the specified observer
%
%            [3,N] = size(pos), double = class(pos)
%
%      lt    the double precision scalar one-way light time or double
%            precision 1xN array of one-way light times between the
%            observer and target in seconds; if the target state
%            is corrected for aberrations, then 'lt' is the
%            one-way light time between the observer and the light
%            time corrected target location
%
%            [1,N] = size(lt), double = class(lt)
%
%            'pos' and 'lt' return with the same vectorization
%            measure (N) as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      %  Load a set of kernels: an SPK file, a PCK
%      %  file and a leapseconds file. Use a meta
%      %  kernel for convenience.
%      %
%      cspice_furnsh( 'standard.tm' )
%
%      %
%      %  Define parameters for a position lookup:
%      %
%      %  Return the position vector of Mars (499) as seen from
%      %  Earth (399) in the J2000 frame
%      %  using aberration correction LT+S (light time plus
%      %  stellar aberration) at the epoch
%      %  July 4, 2003 11:00 AM PST.
%      %
%      target   = 'Mars';
%      epoch    = 'July 4, 2003 11:00 AM PST';
%      frame    = 'J2000';
%      abcorr   = 'LT+S';
%      observer = 'Earth';
%
%      %
%      %  Convert the epoch to ephemeris time.
%      %
%      et = cspice_str2et( epoch );
%
%      %
%      %  Look-up the position for the defined parameters.
%      %
%      [ pos, ltime ] = cspice_spkpos( target, et, frame, ...
%                                        abcorr, observer);
%
%      %
%      %  Output...
%      %
%      txt = sprintf( 'The position of    : %s', target);
%      disp( txt )
%
%      txt = sprintf( 'As observed from   : %s', observer );
%      disp( txt )
%
%      txt = sprintf( 'In reference frame : %s', frame );
%      disp( txt )
%      disp( ' ' )
%
%      txt = sprintf( 'Scalar' );
%      disp( txt )
%
%      utc_epoch = cspice_et2utc( et, 'C', 3 );
%
%      txt = sprintf(  'At epoch           : %s', epoch );
%      disp( txt )
%
%      txt = sprintf(  '                   : i.e. %s', utc_epoch );
%      disp( txt )
%
%      txt = sprintf( ['R (kilometers)     : ' ...
%                      '%12.4f %12.4f %12.4f'], pos );
%      disp( txt )
%
%      txt = sprintf( 'Light time (secs)  : %12.7f', ltime );
%      disp( txt )
%
%      disp(' between observer' )
%      disp(' and target' )
%      disp( ' ' )
%
%      %
%      % Create a vector of et's, starting at 'epoch'
%      % in steps of 100000 ephemeris seconds.
%      %
%      vec_et = [0:4]*100000. + et;
%
%      disp( 'Vector' )
%      vec_epoch = cspice_et2utc( vec_et, 'C', 3 );
%
%      %
%      % Look up the position vectors and light time values
%      % 'ltime'  corresponding to the vector of input
%      % ephemeris time 'vec_et'.
%      %
%      [pos , ltime] = cspice_spkpos( target, vec_et, ...
%                                       frame, abcorr, observer );
%
%      for i=1:5
%
%         txt = sprintf(  'At epoch (UTC)     : %s', vec_epoch(i,:) );
%         disp( txt )
%
%         txt = sprintf( ['R (kilometers)     : ' ...
%                         '%12.4f %12.4f %12.4f'], pos(i) );
%         disp( txt )
%
%         txt = sprintf( ['Light time (secs)  : ' ...
%                        '%12.7f'], ltime(i) );
%         disp( txt )
%
%         disp(' between observer' )
%         disp(' and target' )
%         disp( ' ' )
%
%      end
%
%      %
%      %  It's always good form to unload kernels after use,
%      %  particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      The position of    : Mars
%      As observed from   : Earth
%      In reference frame : J2000
%
%      Scalar
%      At epoch           : July 4, 2003 11:00 AM PST
%                         : i.e. 2003 JUL 04 19:00:00.000
%      R (kilometers)     : 73822235.3105 -27127918.9985 -18741306.3015
%      Light time (secs)  :  269.6898814
%       between observer
%       and target
%
%      Vector
%      At epoch (UTC)     : 2003 JUL 04 19:00:00.000
%      R (kilometers)     : 73822235.3105 -27127918.9985 -18741306.3015
%      Light time (secs)  :  269.6898814
%       between observer
%       and target
%
%      At epoch (UTC)     : 2003 JUL 05 22:46:40.000
%      R (kilometers)     : 73140185.4144 -26390524.7797 -18446763.0348
%      Light time (secs)  :  266.5640394
%       between observer
%       and target
%
%      At epoch (UTC)     : 2003 JUL 07 02:33:20.000
%      R (kilometers)     : 72456239.6608 -25681031.0146 -18163339.1448
%      Light time (secs)  :  263.4803533
%       between observer
%       and target
%
%      At epoch (UTC)     : 2003 JUL 08 06:20:00.000
%      R (kilometers)     : 71771127.0087 -24999259.4606 -17890946.6362
%      Light time (secs)  :  260.4395234
%       between observer
%       and target
%
%      At epoch (UTC)     : 2003 JUL 09 10:06:40.000
%      R (kilometers)     : 71085543.8280 -24345021.1811 -17629490.7100
%      Light time (secs)  :  257.4422002
%       between observer
%       and target
%
%-Particulars
%
%   A sister version of this routine exists named mice_spkpos that returns
%   the output arguments as fields in a single structure.
%
%   Aberration corrections 
%   ====================== 
% 
%   In space science or engineering applications one frequently 
%   wishes to know where to point a remote sensing instrument, such 
%   as an optical camera or radio antenna, in order to observe or 
%   otherwise receive radiation from a target. This pointing problem 
%   is complicated by the finite speed of light:  one needs to point 
%   to where the target appears to be as opposed to where it actually 
%   is at the epoch of observation. We use the adjectives 
%   "geometric," "uncorrected," or "true" to refer to an actual 
%   position or state of a target at a specified epoch. When a 
%   geometric position or state vector is modified to reflect how it 
%   appears to an observer, we describe that vector by any of the 
%   terms "apparent," "corrected," "aberration corrected," or "light 
%   time and stellar aberration corrected." The SPICE Toolkit can
%   correct for two phenomena affecting the apparent location of an
%   object:  one-way light time (also called "planetary aberration") and
%   stellar aberration.
%
%   One-way light time
%   ------------------
%
%   Correcting for one-way light time is done by computing, given an
%   observer and observation epoch, where a target was when the observed
%   photons departed the target's location. The vector from the
%   observer to this computed target location is called a "light time
%   corrected" vector. The light time correction depends on the motion
%   of the target relative to the solar system barycenter, but it is
%   independent of the velocity of the observer relative to the solar
%   system barycenter. Relativistic effects such as light bending and
%   gravitational delay are not accounted for in the light time
%   correction performed by this routine.
% 
%   Stellar aberration
%   ------------------
%
%   The velocity of the observer also affects the apparent location 
%   of a target:  photons arriving at the observer are subject to a 
%   "raindrop effect" whereby their velocity relative to the observer 
%   is, using a Newtonian approximation, the photons' velocity 
%   relative to the solar system barycenter minus the velocity of the 
%   observer relative to the solar system barycenter. This effect is 
%   called "stellar aberration."  Stellar aberration is independent 
%   of the velocity of the target. The stellar aberration formula 
%   used by this routine does not include (the much smaller)
%   relativistic effects.
% 
%   Stellar aberration corrections are applied after light time 
%   corrections:  the light time corrected target position vector is  
%   used as an input to the stellar aberration correction. 
% 
%   When light time and stellar aberration corrections are both 
%   applied to a geometric position vector, the resulting position  
%   vector indicates where the target "appears to be" from the 
%   observer's location.  
% 
%   As opposed to computing the apparent position of a target, one 
%   may wish to compute the pointing direction required for 
%   transmission of photons to the target. This also requires correction 
%   of the geometric target position for the effects of light time 
%   and stellar aberration, but in this case the corrections are 
%   computed for radiation traveling *from* the observer to the target. 
%   We will refer to this situation as the "transmission" case.
%
%   The "transmission" light time correction yields the target's 
%   location as it will be when photons emitted from the observer's 
%   location at `et' arrive at the target. The transmission stellar 
%   aberration correction is the inverse of the traditional stellar 
%   aberration correction:  it indicates the direction in which 
%   radiation should be emitted so that, using a Newtonian 
%   approximation, the sum of the velocity of the radiation relative 
%   to the observer and of the observer's velocity, relative to the  
%   solar system barycenter, yields a velocity vector that points in  
%   the direction of the light time corrected position of the target. 
%  
%   One may object to using the term "observer" in the transmission
%   case, in which radiation is emitted from the observer's location.
%   The terminology was retained for consistency with earlier
%   documentation.
% 
%   Below, we indicate the aberration corrections to use for some 
%   common applications: 
% 
%      1) Find the apparent direction of a target for a remote-sensing
%         observation.
%
%            Use 'LT+S' or 'CN+S: apply both light time and stellar 
%            aberration corrections.
%
%         Note that using light time corrections alone ('LT' or 'CN')
%         is generally not a good way to obtain an approximation to
%         an apparent target vector: since light time and stellar
%         aberration corrections often partially cancel each other,
%         it may be more accurate to use no correction at all than to
%         use light time alone.
%
%
%      2) Find the corrected pointing direction to radiate a signal
%         to a target. This computation is often applicable for 
%         implementing communications sessions.
%
%            Use 'XLT+S' or 'XCN+S: apply both light time and stellar 
%            aberration corrections for transmission.
%
%
%      3) Compute the apparent position of a target body relative
%         to a star or other distant object.
%
%            Use 'LT', 'CN', 'LT+S', or 'CN+S' as needed to match the
%            correction applied to the position of the distant
%            object. For example, if a star position is obtained from
%            a catalog, the position vector may not be corrected for
%            stellar aberration. In this case, to find the angular
%            separation of the star and the limb of a planet, the
%            vector from the observer to the planet should be
%            corrected for light time but not stellar aberration.
%
%
%      4) Obtain an uncorrected state vector derived directly from 
%         data in an SPK file.
%
%            Use 'NONE'.
%
%
%      5) Use a geometric state vector as a low-accuracy estimate
%         of the apparent state for an application where execution 
%         speed is critical.
%
%            Use 'NONE'.
%
%
%      6) While this routine cannot perform the relativistic
%         aberration corrections required to compute states
%         with the highest possible accuracy, it can supply the
%         geometric states required as inputs to these computations.
%
%            Use 'NONE', then apply relativistic aberration
%            corrections (not available in the SPICE Toolkit).
%
% 
%   Below, we discuss in more detail how the aberration corrections 
%   applied by this routine are computed.      
% 
%      Geometric case 
%      ============== 
% 
%      spkezr_c begins by computing the geometric position T(et) of the 
%      target body relative to the solar system barycenter (SSB). 
%      Subtracting the geometric position of the observer O(et) gives 
%      the geometric position of the target body relative to the 
%      observer. The one-way light time, lt, is given by 
% 
%                | T(et) - O(et) | 
%         lt = ------------------- 
%                        c 
% 
%      The geometric relationship between the observer, target, and 
%      solar system barycenter is as shown: 
% 
% 
%         SSB ---> O(et) 
%          |      / 
%          |     / 
%          |    /                            
%          |   /  T(et) - O(et)   
%          V  V                                   
%         T(et) 
% 
% 
%      The returned state consists of the position vector 
% 
%         T(et) - O(et) 
% 
%      and a velocity obtained by taking the difference of the 
%      corresponding velocities. In the geometric case, the  
%      returned velocity is actually the time derivative of the  
%      position. 
% 
% 
%      Reception case 
%      ============== 
% 
%      When any of the options "LT", "CN", "LT+S", "CN+S" is selected
%      for `abcorr', spkezr_c computes the position of the target body at
%      epoch et-lt, where `lt' is the one-way light time. Let T(t) and
%      O(t) represent the positions of the target and observer 
%      relative to the solar system barycenter at time t; then `lt' is 
%      the solution of the light-time equation 
% 
%                | T(et-lt) - O(et) | 
%         lt = ------------------------                            (1) 
%                         c 
% 
%      The ratio  
% 
%          | T(et) - O(et) | 
%        ---------------------                                     (2) 
%                  c 
% 
%      is used as a first approximation to `lt'; inserting (2) into the
%      right hand side of the light-time equation (1) yields the
%      "one-iteration" estimate of the one-way light time ("LT").
%      Repeating the process until the estimates of `lt' converge yields
%      the "converged Newtonian" light time estimate ("CN").
%      
%      Subtracting the geometric position of the observer O(et) gives 
%      the position of the target body relative to the observer: 
%      T(et-lt) - O(et). 
% 
%         SSB ---> O(et) 
%          | \     | 
%          |  \    | 
%          |   \   | T(et-lt) - O(et) 
%          |    \  | 
%          V     V V 
%         T(et)  T(et-lt) 
%       
%      The position component of the light time corrected state  
%      is the vector 
% 
%         T(et-lt) - O(et) 
% 
%      The velocity component of the light time corrected state  
%      is the difference 
% 
%         T_vel(et-lt)*(1-d(lt)/d(et)) - O_vel(et) 
%           
%      where T_vel and O_vel are, respectively, the velocities of the 
%      target and observer relative to the solar system barycenter at 
%      the epochs et-lt and `et'.  
%
%      If correction for stellar aberration is requested, the target 
%      position is rotated toward the solar system 
%      barycenter-relative velocity vector of the observer. The 
%      rotation is computed as follows: 
% 
%         Let r be the light time corrected vector from the observer 
%         to the object, and v be the velocity of the observer with 
%         respect to the solar system barycenter. Let w be the angle 
%         between them. The aberration angle phi is given by 
% 
%            sin(phi) = v sin(w) / c 
% 
%         Let h be the vector given by the cross product 
% 
%            h = r X v 
% 
%         Rotate r by phi radians about h to obtain the apparent 
%         position of the object. 
% 
%      When stellar aberration corrections are used, the rate of change
%      of the stellar aberration correction is accounted for in the
%      computation of the output velocity.
% 
% 
%      Transmission case 
%      ================== 
% 
%      When any of the options "XLT", "XCN", "XLT+S", "XCN+S" is 
%      selected, spkezr_c computes the position of the target body T at 
%      epoch et+lt, where `lt' is the one-way light time. `lt' is the 
%      solution of the light-time equation 
% 
%                | T(et+lt) - O(et) | 
%         lt = ------------------------                            (3) 
%                          c 
% 
%      Subtracting the geometric position of the observer, O(et), 
%      gives the position of the target body relative to the 
%      observer: T(et-lt) - O(et). 
% 
%                 SSB --> O(et) 
%                / |    *  
%               /  |  *  T(et+lt) - O(et)   
%              /   |*      
%             /   *|     
%            V  V  V      
%        T(et+lt)  T(et)     
% 
%      The position component of the light-time corrected state  
%      is the vector 
% 
%         T(et+lt) - O(et) 
% 
%      The velocity component of the light-time corrected state  
%      consists of the difference 
% 
%         T_vel(et+lt)*(1+d(lt)/d(et)) - O_vel(et) 
%           
%      where T_vel and O_vel are, respectively, the velocities of the 
%      target and observer relative to the solar system barycenter at 
%      the epochs et+lt and `et'. 
% 
%      If correction for stellar aberration is requested, the target 
%      position is rotated away from the solar system barycenter- 
%      relative velocity vector of the observer. The rotation is 
%      computed as in the reception case, but the sign of the 
%      rotation angle is negated. 
% 
% 
%   Precision of light time corrections 
%   =================================== 
% 
%      Corrections using one iteration of the light time solution 
%      ---------------------------------------------------------- 
% 
%      When the requested aberration correction is "LT", "LT+S", 
%      "XLT", or "XLT+S", only one iteration is performed in the 
%      algorithm used to compute lt. 
% 
%      The relative error in this computation 
% 
%         | LT_ACTUAL - LT_COMPUTED |  /  LT_ACTUAL 
% 
%      is at most  
% 
%          (V/C)**2 
%         ---------- 
%          1 - (V/C) 
% 
%      which is well approximated by (V/C)**2, where V is the 
%      velocity of the target relative to an inertial frame and C is 
%      the speed of light. 
% 
%      For nearly all objects in the solar system V is less than 60
%      km/sec. The value of C is ~300000 km/sec. Thus the
%      one-iteration solution for LT has a potential relative error
%      of not more than 4e-8. This is a potential light time error of
%      approximately 2e-5 seconds per astronomical unit of distance
%      separating the observer and target. Given the bound on V cited
%      above:
%
%         As long as the observer and target are separated by less
%         than 50 astronomical units, the error in the light time
%         returned using the one-iteration light time corrections is
%         less than 1 millisecond.
%
%         The magnitude of the corresponding position error, given
%         the above assumptions, may be as large as (V/C)**2 * the
%         distance between the observer and the uncorrected target
%         position: 300 km or equivalently 6 km/AU.
%
%      In practice, the difference between positions obtained using
%      one-iteration and converged light time is usually much smaller
%      than the value computed above and can be insignificant. For
%      example, for the spacecraft Mars Reconnaissance Orbiter and
%      Mars Express, the position error for the one-iteration light
%      time correction, applied to the spacecraft-to-Mars center
%      vector, is at the 1 cm level.
%      
%      Comparison of results obtained using the one-iteration and
%      converged light time solutions is recommended when adequacy of
%      the one-iteration solution is in doubt.
%
%
%      Converged corrections 
%      ---------------------
%
%      When the requested aberration correction is 'CN', 'CN+S',
%      'XCN', or 'XCN+S', as many iterations as are required for
%      convergence are performed in the computation of LT. Usually
%      the solution is found after three iterations. The relative
%      error present in this case is at most
%
%          (V/C)**4
%         ----------
%          1 - (V/C)
%
%      which is well approximated by (V/C)**4. 
%
%         The precision of this computation (ignoring round-off
%         error) is better than 4e-11 seconds for any pair of objects
%         less than 50 AU apart, and having speed relative to the
%         solar system barycenter less than 60 km/s.
%
%         The magnitude of the corresponding position error, given
%         the above assumptions, may be as large as (V/C)**4 * the
%         distance between the observer and the uncorrected target
%         position: 1.2 cm at 50 AU or equivalently 0.24 mm/AU.
%
%      However, to very accurately model the light time between
%      target and observer one must take into account effects due to
%      general relativity. These may be as high as a few hundredths
%      of a millisecond for some objects.
%
% 
%   Relativistic Corrections 
%   ========================= 
% 
%   This routine does not attempt to perform either general or 
%   special relativistic corrections in computing the various 
%   aberration corrections. For many applications relativistic 
%   corrections are not worth the expense of added computation 
%   cycles. If however, your application requires these additional 
%   corrections we suggest you consult the astronomical almanac (page 
%   B36) for a discussion of how to carry out these corrections. 
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine spkpos_c.
%
%   MICE.REQ
%   SPK.REQ
%   NAIF_IDS.REQ
%   FRAMES.REQ
%   TIME.REQ
%
%-Version
%
%-  -Mice Version 1.0.2, 03-JUL-2014 (NJB) (BVS) (EDW)
%
%      Discussion of light time corrections was updated. Assertions
%      that converged light time corrections are unlikely to be
%      useful were removed.
%
%   Last update was 07-NOV-2013 (EDW)
%
%      Added aberration algorithm explanation to Particulars section.
%
%   -Mice Version 1.0.1, 22-DEC-2008, EDW (JPL)
%
%      Header edits performed to improve argument descriptions.
%      These descriptions should now closely match the descriptions
%      in the corresponding CSPICE routine.
%
%      Corrected typo in I/O section. Replaced the "ptarg"
%      return argument name with "pos."
%
%   -Mice Version 1.0.0, 22-NOV-2005, EDW (JPL)
%
%-Index_Entries
%
%   using names get target position relative to an observer
%   position relative to observer corrected for aberrations
%   read ephemeris data
%   read trajectory data
%
%-&

function [pos, lt] = cspice_spkpos(targ, et, ref, abcorr, obs)

   switch nargin
      case 5

         targ   = zzmice_str(targ);
         et     = zzmice_dp(et);
         ref    = zzmice_str(ref);
         abcorr = zzmice_str(abcorr);
         obs    = zzmice_str(obs);

      otherwise

         error ( ['Usage: [_pos(3)_, _lt_] = ' ...
                  'cspice_spkpos( `targ`, _et_, `ref`, `abcorr`, `obs`)'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [ptarg] = mice('spkpos_s', targ, et,ref, abcorr, obs);
      pos     = reshape( [ptarg.pos], 3, [] );
      lt      = reshape( [ptarg.lt],  1, [] );
   catch
      rethrow(lasterror)
   end




