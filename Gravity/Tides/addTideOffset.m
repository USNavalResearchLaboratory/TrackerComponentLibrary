function [CTide,STide]=addTideOffset(C,S,TT1,TT2,deltaTTUT1,xpyp,dXdY)
%%ADDTIDEOFFSET Given a set of fully normalized spherical harmonic
%               gravitational coefficients for a tide-free model of the
%               Earth, compute the coefficients for a model with tides
%               including oceanic, solid Earth and pole tides.
%
%INPUTS: C An array holding the fully normalized coefficient terms that are
%          multiplied by cosines in the spherical harmonic expansion. If
%          passed to a CountingClusterSet class, C(n+1,m+1) is the
%          coefficient of degree n and order m. When a maximum degree of M
%          is used, all C must have values for all n from 0 to M and for
%          all m from 0 to n for each n. If coefficients are not present
%          for certain degrees, then insert a 0.
%        S An array class holding the coefficient terms that are multiplied
%          by sines in the harmonic expansion. The requirements on S are
%          the same as those on C.
%  TT1,TT2 Two parts of a Julian date given in terrestrial time (TT). The
%          units of the date are days. The full date is the sum of both
%          terms. The date is broken into two parts to provide more bits of
%          precision. It does not matter how the date is split.
% deltaTTUT1 An optional parameter specifying the difference between TT and
%          UT1 in seconds. Values of this parameter can be obtained from
%          http://www.iers.org/nn_11474/IERS/EN/DataProducts/EarthOrientationData/eop.html?__nnn=true
%          or 
%          http://www.usno.navy.mil/USNO/earth-orientation/eo-products
%          If this parameter is omitted or an empty matrix is passed, then
%          the value provided by the function getEOP will be used instead.
%     xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
%          including the effects of tides and librations. If this parameter
%          is omitted or if an empty matrix is passed, the value from the
%          function getEOP will be used.
%     dXdY dXdY=[dX;dY] are the celestial pole offsets with respect to the
%          IAU 2006/2000A precession/nutation model in radians If this
%          parameter is omitted or an empty matrix is passed, the value
%          from the function getEOP will be used.
%
%OUTPUTS: CTide,STide The same as C and S but modified with the effects of
%                     the oceanic tides.
%
%This function combines the offsets for the coefficients due to tides as
%described in [1]. The function readJPLEphem is used to determine the
%locations of the Sun and the Moon in ICRS coordinates as needed for the
%function gravSolidTideOffset, where the function TT2TDB provides the
%approximate conversion to TDB for readJPLEphem.
%
%To determine the precision of the method, the comments for the functions 
%gravPoleTideOffset, gravSolidTideOffset, and gravOceanTideOffset should be
%consulted.
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5)
    deltaTTUT1=[];
end

if(nargin<6)
    xpyp=[];
end

if(nargin<7)
    dXdY=[];
end

CTide=C;
STide=S;

numEl=length(CTide);

%Compute the offsets due to solid Earth tides.

[TDB1,TDB2]=TT2TDB(TT1,TT2);%Get approximate TDB.
%Sun position with respect to Earth.
SunGCRSPosVel=readJPLEphem(TDB1,TDB2,11,3);
rSunITRS=GCRS2ITRS(SunGCRSPosVel(1:3),TT1,TT2,deltaTTUT1,xpyp,dXdY);

%Moon position with respect to Earth.
MoonGCRSPosVel=readJPLEphem(TDB1,TDB2,10,3);
rMoonITRS=GCRS2ITRS(MoonGCRSPosVel(1:3),TT1,TT2,deltaTTUT1,xpyp,dXdY);

[deltaC,deltaS]=gravSolidTideOffset(rMoonITRS,rSunITRS,TT1,TT2);

%Add in the solid Earth tide offsets.
numOffset=length(deltaC);
if(numOffset<=numEl)
    CTide(1:numOffset)=CTide(1:numOffset)+deltaC(:);
    STide(1:numOffset)=STide(1:numOffset)+deltaS(:);
else
    CTide(:)=CTide(:)+deltaC(1:numEl);
    STide(:)=STide(:)+deltaS(1:numEl);
end

%Compute the offsets due to oceanic tides
[deltaC,deltaS]=gravOceanTideOffset(TT1,TT2);

%Add in the oceanic tide offsets.
numOffset=length(deltaC);
if(numOffset<=numEl)
    CTide(1:numOffset)=CTide(1:numOffset)+deltaC(:);
    STide(1:numOffset)=STide(1:numOffset)+deltaS(:);
else
    CTide(:)=CTide(:)+deltaC(1:numEl);
    STide(:)=STide(:)+deltaS(1:numEl);
end

%Compute the pole tide offsets.
[deltaC,deltaS]=gravPoleTideOffset(TT1,TT2,xpyp);

%Add in the pole tide offsets.
numOffset=length(deltaC.clusterEls);
if(numOffset<=numEl)
    CTide(1:numOffset)=CTide(1:numOffset)+deltaC(:);
    STide(1:numOffset)=STide(1:numOffset)+deltaS(:);
else
    CTide(:)=CTide(:)+deltaC(1:numEl);
    STide(:)=STide(:)+deltaS(1:numEl);
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
