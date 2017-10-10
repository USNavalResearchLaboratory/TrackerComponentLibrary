function [deltaC,deltaS]=gravPoleTideOffset(TT1,TT2,xpyp)
%%GRAVPOLETIDEOFFSET Compute the offsets to add the fully normalized
%                    spherical harmonic gravitational coefficients to
%                    handle the effects of the pole tides. Additional Earth
%                    tides and other effects need to be added separately to
%                    get the full offset due to tides.
%
%%INPUTS TT1,TT2 Two parts of a Julian date given in terrestrial time (TT). 
%                The units of the date are days. The full date is the sum
%                of both terms. The date is broken into two parts to
%                provide more bits of precision. It does not matter how the
%                date is split.
%           xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
%                including the effects of tides and librations. If this
%                parameter is omitted or if an empty matrix is passed, the
%                value from the function getEOP will be used.
%
%OUTPUTS:deltaC An array holding the offsets to the fully normalized
%               coefficient terms that are multiplied by cosines in the
%               spherical harmonic expansion of the gravitational
%               potential. If given to a CountingClusterSet class, then
%               C(n+1,m+1) is the coefficient of degree n and order m. The
%               offsets only go up to the degree and order of the FES2004
%               coefficients.
%        deltaS An array holding the coefficient terms that are multiplied
%               by sines in the harmonic expansion. The format of S is the
%               same as that of C.
%
%This implements the pole tide models that are described in Sections 6.4
%and 6.5 of [1]. When implementing the oceanic pole tide model, only the
%most significant terms are used (the ones highlighted in the conventions)
%rather than the full set of parameters from Desai.
%
%To get the full offset in the gravitational coefficients due to tides,
%the functions gravSolidTideOffset, gravOceanTideOffset and
%gravPoleTideOffset have to be combined.
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5)
    [JulUTC1,JulUTC2]=TT2UTC(TT1,TT2);
    xpyp=getEOP(JulUTC1,JulUTC2);
end

xp=xpyp(1);
yp=xpyp(2);

M=2;
totalNumCoeffs=(M+1)*(M+2)/2;
emptyData=zeros(totalNumCoeffs,1);
deltaC=emptyData;
deltaS=emptyData;

xpBarypBar=meanRotPoleLoc(TT1,TT2);
xpBar=xpBarypBar(1);
ypBar=xpBarypBar(2);

%From Equation 7.24 in Section 7.1.4
m1=xp-xpBar;
m2=-(yp-ypBar);

rad2ArcSec=(180/pi)*60*60;
m1=m1*rad2ArcSec;
m2=m2*rad2ArcSec;

%The solid Earth pole tide from Section 6.4
%deltaC_{2,1}
deltaC(4)=-1.333e-9*(m1+0.115*m2);
%deltaS_{2,1}
deltaS(4)=-1.333e-9*(m2-0.115*m1);

%The main coefficients from the ocean pole tide in Section 6.5
%deltaC_{2,1}
deltaC(4)=deltaC(4)+-2.1778e-10*(m1-0.01724*m2);
%deltaS_{2,1}
deltaS(4)=deltaS(4)+-1.7232e-10*(m2-0.03365*m1);
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
