function [accel,C,S]=EGM08EarthAccel(rVec,accelOptions,M,TT1,TT2,effectsToInclude,EOP,C,S)
%%EGM08EARTHACCEL Get acceleration due to gravity from the Earth. This can
%            be a simple J2 model in a generic (not rigorously defined)
%            Earth-centered Earth-fixed (ECEF) or Earth-centered inertial
%            (ECI) coordinate system, or it can be more rigorously
%            defined in the International Terrestrial Reference System
%            (ITRS) or the Geocentric Celestial Reference System (GCRS),
%            accounting for effects such as polar motion, tides, and the
%            drift of the coefficients over time. The EGM2008 acceleration
%            model is used for parameters.
%
%INPUTS: rVec A 3XN (or 6XN) set of vectors of position (and velocity) in
%             meters in the same coordinate system as desired for the
%             output as specified by the parameter accelOptions. This is
%             the set of N positions where the acceleration due to the
%             Earth's gravity is desired. The velocity component is only
%             used if accelOptions specified an Earth-fixed coordinate
%             system. Specifically, if accelOptions is 0, 3, or 5. If a 3XN
%             vector is provided, but the velocity component is required,
%             the extra three elements are taken to be all zeros.
% accelOptions An optional parameter indicating the coordinate system and
%             options for the acceleration. Possible values are
%             0 (The default if omitted) The acceleration is found in a
%               generic ECEF coordinate system using either a Keplerian
%               model (if M=1) the J2 gravitational model (if M=2, the
%               default) with J2 given by the corresponding coefficient
%               from the zero-tide EGM2008 model without correcting for
%               polar motion, other tides, or the drift of the coefficients
%               over time. It is assumed that the rotation axis is aligned
%               with the z axis. Being an accelerating coordinate system,
%               Newtonian corrections for the Coriolis effect are used.
%               Values M>=3 are invalid. No Earth orientation parameters 
%               are used. The time is also not used and neither is the
%               effectsToInclude nor C and S.
%             1 The same as 0 but in a generic Earth-centered inertial
%               coordinate system, meaning that the Coriolis effect is not
%               taken into account.
%             2 The acceleration is found using the first M elements of the
%               EGM2008 model in GCRS coordinates. The additional effects
%               accounted for (beyond just using a zero-tide model) are
%               determined by the effectsToInclude boolean vector. If C and
%               S are not provided, then the getEGMGravCoeffs function is
%               used. Earth orientation parameters are used, including the
%               length-of-day (LOD) parameter for the rotation rate of the
%               Earth.
%             3 The same as 2 but in ITRS coordinates, meaning that a
%               Newtonian Coriolis term is added to the coefficients,
%               because it is an accelerating coordinate system.
%           M The number of terms in the EGM2008 model to include. If
%             accelOptions=0 or accelOptions=1, then M chooses between
%             Keplerian, or J2 models. Otherwise, the number indicates
%             the highest order of the full EGM2008 model used. If omitted,
%             or an empty matrix is passed, a default value of 2 is used
%             for accelOptions=0 or 1 and a default of 3 is used for
%             accelOptions=2 or 3. If Inf or another number larger than the
%             highest coefficient order in the EGM2008 model is passed,
%             then the total number of coefficients in the EGM2008 model is
%             used.
%     TT1,TT2 Two parts of a Julian date given in terrestrial time (TT).
%             The units of the date are days. The full date is the sum of
%             both terms. The date is broken into two parts to provide more
%             bits of precision. It does not matter how the date is split.
%             The date is only used for accelOptions=2-3. Otherwise it can
%             be omitted or empty matrices can be passed.
% effectsToInclude A boolean vector indicating which effects are to be
%             included (if accelOptions>1). If an element of the vector is
%             1, that means that the effect will be taken into account. If
%             omitted, the default is [1;1;0;0;0]; The elements of the
%             vector are:
%             1) The drift of the low-order coefficients over time, as
%                given by the getGravCoeffOffset4Drift function.
%             2) The effects of polar motion on the coefficients, as given
%                by the getAdjustedGravCoeffs4PolarMotion function. They
%                are added after tidal effects (if any).
%             3) The effects of solid Earth tides as given by the 
%                gravSolidTideOffset function.
%             4) The effects of pole tides as given by the
%                gravPoleTideOffset function.
%             5) The effects of ocean tides as given by the
%                gravOceanTideOffset function.
%         EOP A structure containing the Earth orientation parameters at
%             the given time. If omitted or an empty matrix is passed, the
%             values from the getEOP function are used. These are not used
%             if accelOptions=0 or 1. The elements of the structure are:
%             xpyp These are the polar motion coordinates in radians
%                  including the effects of tides and librations.
%             dXdY The celestial pole offsets with respect to the IAU
%                  2006/2000A precession/nutation model in radians.
%             deltaTTUT1 The difference between TT and UT1 in seconds
%              LOD The difference between the length of the day using TT
%                  and the length of the day in UT1. This is an
%                  instantaneous parameter (a derivative). The units are
%                  seconds.
%              C,S Optionally, the TIDE FREE EGM2008 coefficients with
%                  adjType=0 (the default) as obtained from the
%                  getEGMGravCoeffs function. Passing these saves a call to
%                  the getEGMGravCoeffs function and speeds up multiple
%                  calls to this function. These parameters are ignored if
%                  accelOptions=0 or 1.
%
%OUTPUTS: accel A 3XN matrix of the N accelerations due to gravity of the
%               Earth in meters per second squared in the specified
%               coordinate system (one for each input rVec).
%           C,S Values of the tide-free EGM2008 spherical harmonic
%               coefficients that can be passed to subsequent calls of this
%               function with the same options). Note that C and S are
%               modified (and then changed back), so care must be taken
%               with multithreading. Emptyt matrices are returned if
%               accelOptions=0 or accelOptions=1;
%
%For simple models, many parameters can be omitted.
%The J2 model in an ECEF system can be obtained using
%accel=EGM08EarthAccel(rVec)
%The function can be used for the Keplerian model in an ECEF system using
%accel=EGM08EarthAccel(rVec,0,1)
%or in an ECI system using
%accel=EGM08EarthAccel(rVec,1,1)
%The J2 model can be obtained in an ECI system using
%accel=EGM08EarthAccel(rVec,1,2)
%
%The J2 gravitational model is derived in Appendix E of [1].
%The Coriolis and centrifugal force models due to the use of non-inertial
%coordinate systems are derived from Appendix A of the same document.
%The various other effects are discussed in other sections.
%
%REFERENCES:
%[1] D. F. Crouse, "An Overview of Major Terrestrial, Celestial, and
%    Temporal Coordinate Systems for Target Tracking," Formal Report, Naval
%    Research Laboratory, no. NRL/FR/5344--16-10,279, 10 Aug. 2016, 173
%    pages.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%incorporating elements from a Coriolis correction by David Karnick.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(accelOptions))
    accelOptions=0;
end

if(nargin<3||isempty(M))
    if(accelOptions==0||accelOptions==1)
        M=2;
    else
        M=3;
    end
end
%If someone passes a number greater than the total number of coefficients,
%then just limit it to the total number in the EGM2008 model.
M=min(2190,M);

%Universal gravitational constant times the mass of the Earth
GM=Constants.EGM2008GM;

%Semi-major axis of the Earth
a=Constants.EGM2008SemiMajorAxis;

%The rotation rate of the Earth in radians per second.
omega=Constants.EGM2008EarthRotationRate;

%Check for the easy cases Keplerian or J2.
if(accelOptions==0||accelOptions==1)
    r=rVec(1:3,:);%Positions
    if(size(rVec,1)>3)
        v=rVec(4:6,:);%Velocities.
    else
        v=zeros(3,size(rVec,2));
    end
    
    rMag=sqrt(sum(r.*r,1));

    %Newton's law for a point or sphere.
    aNewton=-bsxfun(@times,(GM./rMag.^2),bsxfun(@rdivide,r,rMag));

    switch(M)
        case 1%Simple Keplerian dynamics.
            accel=aNewton;
        case 2%The J2 dynamic model
            %Putting this in explicitly saves a call to getEGMGravCoeffs
            %for just one parameter. This is the zero-tide value.
            C20Bar=-0.484169317366974e-03;
            C20=C20Bar*sqrt(5);
            J2=-C20;
            
            temp=5*r(3,:).^2./rMag.^2;
            %First-order oblateness correction.
            aOblateness=-bsxfun(@times,((3*GM*a.^2*J2./(2*rMag.^5))),[r(1,:).*(1-temp);
                                                                     r(2,:).*(1-temp);
                                                                     r(3,:).*(3-temp)]);     

            accel=aNewton+aOblateness;
        otherwise
            error('Invalid value of M for the chosen accelOptions value')
    end

    %If Coriolis terms should be added due to the rotation of the
    %Earth.
    if(accelOptions==0)
        %The rotation vector for the Earth is
        %Omega=[0;0;1]*omega;
        %and the Coriolis and centrifugal forces to add to accel are
        %aCoriolis=-2*bsxfun(@cross,Omega,v);
        %aCentrifugal=-bsxfun(@cross,Omega,bsxfun(@cross,Omega,r));
        
        %However, it is faster to add the components directly:
        accel(1,:)=accel(1,:)+2*omega*v(2,:)+omega^2*r(1,:);
        accel(2,:)=accel(2,:)-2*omega*v(1,:)+omega^2*r(2,:);
    end
    C=[];
    S=[];
    return;
end

%Get the spherical harmonic coefficients up to the requested order if they
%are not provided.
if(nargin<8||isempty(C))
    isTideFree=true;
    modelType=0;
    [C,S]=getEGMGravCoeffs(M,isTideFree,modelType);
end

if(nargin<6||isempty(effectsToInclude))
    effectsToInclude=[1;1;0;0;0];
end

if(nargin<7||isempty(EOP))
    [JulUTC1,JulUTC2]=TT2UTC(TT1,TT2);
    [xpyp,dXdY,~,deltaTTUT1,LOD]=getEOP(JulUTC1,JulUTC2);
else
    xpyp=EOP.dxdy;
    dXdY=EOP.dXdY;
    deltaTTUT1=EOP.deltaTTUT1;
    LOD=EOP.LOD;
end

%The second one is added after all of the others.
numDelta=sum(effectsToInclude)-effectsToInclude(2);
deltaC=cell(numDelta,1);
deltaS=cell(numDelta,1);

curDelta=1;
if(effectsToInclude(1)~=false)
    %Get the offsets to C and S due to the drift of the coefficients.
    [deltaC{curDelta},deltaS{curDelta}]=getGravCoeffOffset4Drift(TT1,TT2,0);
    curDelta=curDelta+1;
end

if(effectsToInclude(3)~=false)
    %Compute the offsets due to solid Earth tides.
    
    [TDB1,TDB2]=TT2TDB(TT1,TT2);%Get approximate TDB.
    %Sun Position with respect to Earth.
    SunGCRSPosVel=readJPLEphem(TDB1,TDB2,11,3);
    rSunITRS=GCRS2ITRS(SunGCRSPosVel(1:3),TT1,TT2,deltaTTUT1,xpyp,dXdY);

    %Moon position with respect to Earth.
    MoonGCRSPosVel=readJPLEphem(TDB1,TDB2,10,3);
    rMoonITRS=GCRS2ITRS(MoonGCRSPosVel(1:3),TT1,TT2,deltaTTUT1,xpyp,dXdY);

    [deltaC{curDelta},deltaS{curDelta}]=gravSolidTideOffset(rMoonITRS,rSunITRS,TT1,TT2);
    curDelta=curDelta+1;
end

if(effectsToInclude(4)~=false)
    [deltaC{curDelta},deltaS{curDelta}]=gravPoleTideOffset(TT1,TT2,xpyp);
    curDelta=curDelta+1;
end

if(effectsToInclude(5)~=false)
    [deltaC{curDelta},deltaS{curDelta}]=gravOceanTideOffset(TT1,TT2);
    curDelta=curDelta+1;
end

%Find out the total number of coefficients in C and S that need to be
%changed --this is the maximum number of elements in the deltas.
numCoeffChanged=0;
for curDelta=1:numDelta
    curLength=length(deltaC{curDelta});
    if(curLength>numCoeffChanged)
        numCoeffChanged=curLength;
    end
end

if(length(C)<numCoeffChanged)
    numCoeffChanged=length(C);
end

%Now, save the elements in C and S, so that they can be restored for the
%return value when the function exits.
CSaved=C(1:numCoeffChanged);
SSaved=S(1:numCoeffChanged);

%Now, add in all of the effects, except the effects of polar motion on the
%coefficients.
for curDelta=1:numDelta
    numOffset=length(deltaC{curDelta});
    num2Change=min(numOffset,numCoeffChanged);
    
    C(1:num2Change)=C(1:num2Change)+deltaC{curDelta}(1:num2Change);
    S(1:num2Change)=S(1:num2Change)+deltaS{curDelta}(1:num2Change);
end

%Add in the effects of polar motion. This just changes the C21 and S21
%elements, which would otherwise be zero.
if(effectsToInclude(2)~=false)
    [~,~,C,S]=getAdjustedGravCoeffs4PolarMotion(C,S,TT1,TT2,true);
end

%Now, get the acceleration due to gravity in ITRS coordinates WITHOUT the
%Coriolis effect. This can be rotated to GCRS coordinates, if needed or
%have Coriolis terms added, if staying in this coordinate system. This
%requires getting the position in ITRS, spherical coordinates.

if(accelOptions~=3)%If r is not already in ITRS coordinates
   %Convert r (and v if present) into ITRS coordinates 
   rVec=GCRS2ITRS(rVec,TT1,TT2,deltaTTUT1,xpyp,dXdY,LOD);
end

%Convert the position components into spherical coordinates for the
%spherHarmonicEval function.
r=rVec(1:3,:);%Positions
rSpher=Cart2Sphere(r);

[~,accel]=spherHarmonicEval(C,S,rSpher,a,GM);

%Now, if the acceleration is supposed to be in ITRS coordinates, then add
%the Coriolis terms. Otherwise, rotate it into GCRS coordinates.
if(accelOptions==3)
    %The output is in (approximate) ITRS coordinates, compute the Coriolis
    %terms.
    v=rVec(4:6,:);%The velocity in ITRS coordinates.
    
    %Get the axis of rotation in ITRS coordinates. This is the z-axis in
    %TIRS coordinates.
    omegaAxis=TIRS2ITRS([0;0;1],TT2,TT2,xpyp);
    
    %Modify the value of omega to deal with the LOD Earth orientation
    %parameter. This Equation is in Section IIIB of Crouse. The 86400 is
    %the number of seconds in a Julian day.
    omega=omega*(1-LOD/86400);
    
    %The rotation vector for the Earth
    Omega=omegaAxis*omega;

    aCoriolis=-2*bsxfun(@cross,Omega,v);
    aCentrifugal=-bsxfun(@cross,Omega,bsxfun(@cross,Omega,r));
    accel=accel+aCoriolis+aCentrifugal;
else
    %Rotate the output into GCRS coordinates. This is just a rotation (as
    %one would apply to position components). Thus, we do not want to pass
    %it as if it were a velocity input to this function.
    accel=ITRS2GCRS(accel,TT1,TT2,deltaTTUT1,xpyp,dXdY,LOD);
end

%Finally, restore the adjusted elements of C and S in case they are to be
%returned.
C(1:numCoeffChanged)=CSaved;
S(1:numCoeffChanged)=SSaved;
%Undo any changes from the getAdjustedGravCoeffs4PolarMotion function.
C(4)=0;
S(4)=0;

end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
