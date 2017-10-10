function [deltaC,deltaS]=gravSolidTideOffset(rMoon,rSun,TT1,TT2)
%%GRAVSOLIDTIDEOFFSET Compute the offsets to add the fully normalized
%                     spherical harmonic gravitational coefficients to
%                     handle the effects of the solid Earth tides. Ocean
%                     tides and other effects need to be added
%                     separately  to get the full offset due to tides.
%                     These offsets are for fully normalized spherical
%                     harmonic tide-free models, not zero-tide models.
%
%INPTUS: rMoon A vector from the center of the Earth in the International
%          Terrestrial Reference System (the WGS-84 coordinate system), an
%          Earth-centered, Earth-fixed coordinate system, to the
%          moon. This can be obtained using the readJPLEphem and GCRS2ITRS
%          functions.
%     rSun A vector from the center of the Earth in the International
%          Terrestrial Reference System (the WGS-84 coordinate system), an
%          Earth-centered, Earth-fixed coordinate system, to the
%          Sun. This can be obtained using the readJPLEphem and GCRS2ITRS
%          functions.
%  TT1,TT2 Two parts of a Julian date given in terrestrial time (TT). The
%          units of the date are days. The full date is the sum of both
%          terms. The date is broken into two parts to provide more bits of
%          precision. It does not matter how the date is split.
%
%OUTPUTS: deltaC An array holding the offsets to the fully normalized
%                coefficient terms that are multiplied by cosines in the
%                spherical harmonic expansion of the gravitational
%                potential. if given to a CountingClusterSet class, then
%                C(n+1,m+1) is the coefficient of degree n and order m. The
%                offsets only go up to degree and order 4.
%         deltaS An array holding the coefficient terms that are multiplied
%                by sines in the spherical harmonic expansion. The format
%                of S is the same as that of C.
%
%The algorithm is the one described in Chapter 5 of [1], except values for
%k_{2m}^{(+)} taken from the 2010 IERS conventions in [2] for an elastic
%Earth model are added, because [1] does not provide the necessary
%coefficients.
%
%The model provided in [1] is not as high order a model as that used in
%Section 6 of [2]. However, the second step of the model in[2] is not easy
%to understand.
%
%To test this function, one can use the same values for the Sun and Moon
%that are given in the IERS's implementation, DEHANTTIDEINEL.F for a
%particular date. Those are
% rSun=[137859926952.015;54228127881.4350;23509422341.6960];
% rMoon=[-179996231.920342;-312468450.131567;-169288918.592160];
% [Jul1,Jul2]=Cal2UTC(2009,4,13,0,0,0);
% [TT1,TT2]=UTC2TT(Jul1,Jul2);
% [deltaC,deltaS]=gravSolidTideOffset(rMoon,rSun,TT1,TT2);
%
%To get the full offset in the gravitational coefficients due to tides,
%the functions gravSolidTideOffset, gravOceanTideOffset and
%gravPoleTideOffset have to be combined.
%
%REFERENCES:
%[1] S. E. Urban and P. K. Seidelmann, Explanatory Supplement to the
%    Astronomical Almanac, 3rd ed., Eds. Mill Valley, CA: University
%    Science Books, 2013.
%[2] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%The convention for indexing arrays of constant parameters is that 1 refers
%to a value with respect to lunar parameters and 2 refers to a value with
%respect to solar parameters.

%The equatorial radius of the Earth
Re=Constants.EarthEqRadius;
GMRat(1)=Constants.MoonEarthMassRatio;
GMRat(2)=Constants.SunEarthMassRatio;

%The corrections are only up to the fourth order coefficients.
M=4;
totalNumCoeffs=(M+1)*(M+2)/2;
emptyData=zeros(totalNumCoeffs,1);
deltaC=CountingClusterSet(emptyData);
deltaS=CountingClusterSet(emptyData);

%Obtain the spherical coordinates of the moon.
pointSpherMoon=Cart2Sphere(rMoon);
rj(1)=pointSpherMoon(1);
Lambdaj(1)=pointSpherMoon(2);
Phij(1)=pointSpherMoon(3);

%Obtain the spherical coordinates of the sun.
pointSpherSun=Cart2Sphere(rSun);
rj(2)=pointSpherSun(1);
Lambdaj(2)=pointSpherSun(2);
Phij(2)=pointSpherSun(3);

%%%%These are the Step 1 corrections

%The nomial value from Table 5.5
k2=0.3;

%Implementing Equations 5.74, 5.75, and 5.76
sumVal20=0;
sumVal21=0;
sumVal22=0;
for n=1:2
    P=legendre(2,sin(Phij(n)));
    P20=P(1);
    P21=P(2);
    P22=P(3);
    
    sumVal20=sumVal20+GMRat(n)*(Re/rj(n))^3*P20;
    sumVal21=sumVal21+GMRat(n)*(Re/rj(n))^3*P21*exp(-1i*Lambdaj(n));
    sumVal22=sumVal22+GMRat(n)*(Re/rj(n))^3*P22*exp(-1i*2*Lambdaj(n));
end
deltaC(2+1,0+1)=(1/sqrt(5))*k2*sumVal20;

deltaC(2+1,1+1)=(1/3)*sqrt(3/5)*k2*real(sumVal21);
deltaS(2+1,1+1)=-(1/3)*sqrt(3/5)*k2*imag(sumVal21);

deltaC(2+1,2+1)=(1/12)*sqrt(12/5)*k2*real(sumVal22);
deltaS(2+1,2+1)=-(1/12)*sqrt(12/5)*k2*imag(sumVal22);

%In the book, no values for k_{2m}^{(+)} were given, so the values here are
%taken from the elastic Earth model in Table 6.3 of the IERS 2010
%conventions.
knmPlus=[-0.00087;-0.00079;-0.00057];
%This loop implements Equation 5.77
for m=0:2
    
    sumVal=0;
    for n=1:2
        P=legendre(2,sin(Phij(n)));
        P2m=P(m+1);
        
        sumVal=sumVal+GMRat(n)*(Re/rj(n))^3*P2m*exp(-1i*m*Lambdaj(n));
    end
    
    deltaC(4+1,m+1)=(1/5)*knmPlus(m+1)*real(sumVal);
    deltaS(4+1,m+1)=-(1/5)*knmPlus(m+1)*imag(sumVal);
end

%%%%Now for the step 2 (frequency-dependent) corrections.
%Convert the terrestrial time to Greenwich mean sidereal time in radians.
GMST=TT2GMST(TT1,TT2);

%These parameters are taken from Table 5.5. The columns are
%Doodson parameters of tau, s, h, p, N', P1 and Am*delta_kx*Hs
m1TideParam=[1, -1,  0,  0,  0,  0,  -16.4e-12;
             1,  1, -2,  0,  0,  0,  -49.6e-12;
             1,  1,  0,  0, -1,  0,   -9.4e-12;
             1,  1,  0,  0,  0,  0,  507.4e-12;
             1,  1,  0,  0,  1,  0,   73.5e-12;
             1,  1,  1,  0,  0, -1,  -15.2e-12];
m2TideParam=[2,  0,  0,  0,  0,  0,   39.5e-12;
             2,  2,  -2, 0,  0,  0,   18.4e-12];

%Turn terrestrial time into barycentric dynamical time TDB.
[TDB1,TDB2]=TT2TDB(TT1,TT2);

%Get the Delaunay variables in radians
vec=DelaunayVar(TDB1,TDB2);
l=vec(1);
lp=vec(2);
F=vec(3);
D=vec(4);
Omega=vec(5);

s=F+Omega;
beta(2)=s;%s: Moon's mean longitude
beta(1)=GMST+pi-s;%tau: time angle in lunar days reckoned from lower transit
beta(3)=s-D;%h: Sun's mean longitude
beta(4)=s-l;%p: longitude of the Moon's mean perigee
beta(5)=-Omega;%N':negative longitude of the Moon's mean node
beta(6)=s-D-lp;%pl: Longitude of the Sun's mean perigee.

%Diurnal tide terms affect the n=2, m=1 coefficients.
n=2;
m=1;
sumVal=0;
for curThetaS=1:size(m1TideParam,1)
    thetaS=sum(beta.*m1TideParam(curThetaS,1:6));
    
    sumVal=sumVal+m1TideParam(curThetaS,7)*(-1i)*exp(1i*thetaS);
end
deltaC(n+1,m+1)=deltaC(n+1,m+1)+real(sumVal);
deltaS(n+1,m+1)=deltaS(n+1,m+1)-imag(sumVal);

%Semi-durnal tides affects the n=2, m=2 coefficients.
n=2;
m=2;

sumVal=0;
for curThetaS=1:size(m2TideParam,1)
    thetaS=sum(beta.*m2TideParam(curThetaS,1:6));
    
    sumVal=sumVal+m2TideParam(curThetaS,7)*(-1i)*exp(1i*thetaS);
end
deltaC(n+1,m+1)=deltaC(n+1,m+1)+real(sumVal);
deltaS(n+1,m+1)=deltaS(n+1,m+1)-imag(sumVal);

deltaC=deltaC.clusterEls;
deltaS=deltaS.clusterEls;
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
