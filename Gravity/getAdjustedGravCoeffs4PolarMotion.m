function [C21,S21,C,S]=getAdjustedGravCoeffs4PolarMotion(C,S,TT1,TT2,updateCAndS,deltaTTUT1,clockLoc)
%%GETADJUSTEDGRAVCOEFFS4POLARMOTION Get the adjusted spherical harmonic
%                   gravitational coefficients C21 and S21 to account for
%                   polar motion as per the IERS 2010 conventions. Polar
%                   motion is the motion of the motion of the rotation axis
%                   of the Earth with respect to the crust of the Earth.
%                   If the crust moves over the liquid core (the
%                   rotation axis moves, then the distribution of mass with
%                   respect to the crust is a little different. The polar
%                   motion effects should probably be added after the drift
%                   of the gravitation model coefficients has been taken
%                   into account. Based on the description in the IERS
%                   conventions, the ZERO-TIDE value of the C20 term
%                   appears to be the appropriate one to use. The polar
%                   motion should be added before other tides.
%
%INPUTS: C, S Arrays of fully normalized spherical harmonic coefficients
%             for a gravitational potential model like those obtained using
%             the getEGMGravCoeffs function. Only coefficients up to order
%             2 are needed. The coefficients are with respect to a
%             coordinate system fixed on the Earth's crust --the
%             International Terrestrial Reference System (ITRS), which can
%             be used interchangably with the WGS-84 coordinate system. The
%             update uses the C20 parameter, and thus depends on the
%             permanent tide model chosen.
%     TT1,TT2 Two parts of a Julian date given in terrestrial time (TT).
%             The units of the date are days. The full date is the sum of
%             both terms. The date is broken into two parts to provide
%             more bits of precision. It does not matter how the date is
%             split.
% updateCAndS An optional parameter. If true, the appropriate values in C
%             and S will be directly modified with the returned values
%             (which is possible since the ClusterSet is a handle class
%             subset). If false (the default if omitted), C and S are not
%             modified. Note that this function only modified C(2+1,1+1)
%             and S(2+1,1+1) if updateCAndS is true and it does not depend
%             on the current values of those parameters, so this function
%             can be called multiple times to update C and S as time
%             progresses.
%  deltaTTUT1 An optional parameter specifying the offset between TT and
%             UT1 in seconds. If this parameter is omitted or an empty
%             matrix is passed, then the value of the function getEOP will
%             be used.
%    clockLoc An optional 3X1 vector specifying the location of the clock
%             in WGS-84 ECEF Cartesian [x;y;z] coordinates with units of
%             meters. Due to relativistic effects, clocks that are
%             synchronized with respect to TT are not synchronized with
%             respect to TDB. If this parameter is omitted, then a clock at
%             the center of the Earth is used.
%    clockLoc An optional 3X1 vector specifying the location of the clock
%             in WGS-84 ECEF Cartesian [x;y;z] coordinates with units of
%             meters. Due to relativistic effects, clocks that are
%             synchronized with respect to TT are not synchronized with
%             respect to TDB. If this parameter is omitted, then a clock at
%             the center of the Earth is used. The fidelity of the model
%             probably is not high enough for this parameter ot matter.
%
%OUTPUTS: C21, S21 The values that would replace C(2+1,1+1) and S(2+1,1+1)
%             when accounting for polar motion if C and S were placed in a
%             CountingClusterSet class to simplify accessing elements.
%        C, S The values of C and S modified with the adjustments.
%
%The effects of polar motion on gravitational coefficients is decribed on
%pages 80-81  (Chapter 6.1) of [1]. The model does not include Earth-
%orientation parameters for the offset of the pole from the mean model,
%which suggests that the model is not sufficiently high fidelity for that
%to matter.
%
%This function uses the meanRotPoleLoc, which has to call TT2BessEpoch. The
%conversion to a besselian epoch depends on the location of the 
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<7)
    clockLoc=[];
end

if(nargin<6)
    deltaTTUT1=[];
end

if(nargin<5||isempty(updateCAndS))
    updateCAndS=false;
end

xpBarypBar=meanRotPoleLoc(TT1,TT2,deltaTTUT1,clockLoc);
xpBar=xpBarypBar(1);
ypBar=xpBarypBar(2);

C20=C(4);
C22=C(6);
S22=S(6);

%This is Equation 6.5 in the IERS 2010 conventions.
C21= sqrt(3)*xpBar*C20-xpBar*C22+ypBar*S22;
S21=-sqrt(3)*ypBar*C20-ypBar*C22-xpBar*S22;

if(updateCAndS)
    C(5)=C21;
    S(5)=S21;
end

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
