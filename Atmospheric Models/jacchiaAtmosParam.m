function [rho,P,T,Te]=jacchiaAtmosParam(Jul1,Jul2,point,F10,F10b,Kp)
%%JACCHIAATMOSPARAM Get basic parameters for atmospheric density,
%                   pressure and temperature from the Jacchia 1971
%                   atmospheric model. The model is best suited for
%                   altitudes above 90km.
%
%INPUTS: Jul1, Jul2 Two parts of a pseudo-Julian date given in UTC. The
%                   units of the date are days. The full date is the sum of
%                   both terms. The date is broken into two parts to
%                   provide more bits of precision. It does not matter how
%                   the date is split.
%             point The [lat;lon;alt] geodetic location under
%                   consideration with latitude and longitude in radians
%                   and altitude in meters.
%               F10 The actual solar flux at 10.7cm wavelength, with units
%                   of 10^-22 W/(m^2*Hz). This is taken to be the average
%                   over the day before the date under consideration.
%              F10b The average solar flux at 10.7cm wavelength, with units
%                   of 10^-22 W/(m^2*Hz). This is the average over three
%                   solar rotations of 27 days.
%                Kp The three-hourly planetary geomagnetic index for a time
%                   6.7 hours earlier than the time under consideration.
%
%OUTPUTS: rho The atmospheric density at the point in question in units of
%             kilograms per cubic meter.
%           P The atmospheric pressure at the point in question in units of
%             Newtons per square meter (Pascals). It assumes that the
%             gasses can be treated as ideal gasses.
%           T The temperature at the point in question with units of
%             degrees Kelvin.
%          Te The exospheric temperature at the location in question with
%             units of degrees Kelvin.
%
%The 1971 Jacchia atmospheric model computes atmospheric densities in two
%major steps. First, the exospheric temperature is computed based on solar
%and geomagnetic data. Then, the density is calculated using a
%bi-polynomial fit. The algorithms in this function are taken primarily
%from Section 3.5.3 of [1] with guidance and some substitutions from the
%original paper [2]. This function does not take into account the seasonal
%variations in helium or hydrogen concentrations.
%
%The atmospheric pressure is obtained using the Ideal Gas Law as described
%in THE NRLMSISE-00 AND HWM-93 USERS GUIDE: Version 1.50, November 2003 by
%Douglas P. Drob.
%
%The anomalous oxygen parameter in the NRLMSISE-00 model is not used.
%
%The position of the Sun is obtained using the readJPLEphem and GCRS2ITRS
%functions. Default values from getEOP are used for all conversions
%requiring them.
%
%REFERENCES:
%[1] O. Montenbruck and E. Gill, "Satellite orbits: models, methods and
%    applications," Springer Science & Business Media. 2012.
%[2] L. G. Jacchia, "Revised Static Models of the Thermosphere and
%    Exosphere with Empirical Temperature PRofiles," SAO Special Report
%    332, Cambridge, 1971
%
%April 2015, David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%%% Define constants
lat=point(1);
lon=point(2);
Z=point(3)/1000;

cn=[28.82678,-7.40066e-2,-1.19407e-2,4.51103e-4,-8.21895e-6,1.07561e-5,-6.97444e-7];%m-profile coefficients
T0=183;%Boundary condition at 90k
z0=90;%Boundary location
zx=125;%Inflection point

% These constants have potential use with added features from Jacchia
% grav=980.665/((1+Z/6356.766)^2);%acceleration due to gravity
% q0=[0.7811,0.20955,0.0093432,6.1471e-6,1];%Fraction by volume for N2,O2,Ar,He,Sum
% Mi=[28.0134,31.9988,39.948,4.0026];%Molecular weight for N2,O2,Ar,He
% M0=28.960;%Mean molecular mass
% rho0=3.46e-9;%Boundary condition at 90k

%%% Calculate exospheric temperature - constants from J71 model
Tc=379.0 + 3.24*F10b + 1.3*(F10-F10b);
%where F10 is actual solar flux at 10.7cm and F10b is average solar flux at
%the same wavelength.

[TT1,TT2]=UTC2TT(Jul1,Jul2);
[TDB1,TDB2]=TT2TDB(TT1,TT2);

%Calculate position of sun
%Sun position with respect to Earth.
SunGCRSPosVel=readJPLEphem(TDB1,TDB2,11,3);
rITRS=GCRS2ITRS(SunGCRSPosVel(1:3),TT1,TT2,deltaTTUT1,xpyp,dXdY);
rSun=Cart2Ellipse(rITRS);
dec=rSun(1);%declination
rObs=ellips2Cart(point);

LAT=TT2LAT(TT1,TT2,rObs);
lha=LAT-pi;%local hour angle

%Actual temp is dependent on position of sun
beta=-37*pi/180;
p=6*pi/180;
gamma=43*pi/180;
m=2.2;
n=3.0;
R=0.3;

eta = 0.5 * abs(lat - dec);
theta = 0.5 * abs(lat + dec);
tau = lha + beta + p*sin(lha + gamma);
tau=wrapRange(tau,-pi,pi,false);
Tl=Tc*(1 + R*(sin(theta)^m + (cos(eta)^m-sin(theta)^m)*cos(tau/2)^n));

%Component due to geomagnetic activities
Thigh=28*Kp + 0.03*exp(Kp); %alt>350km
Tlow=14*Kp + 0.02*exp(Kp); %alt<350km

%A transition function to retain continuity between high and low alt
f=0.5*(tanh(0.04*(Z-350))+1);
Tg=f*Thigh+(1-f)*Tlow;

%Note: Jacchia '71 determined that an additional temperature component due 
%to semiannual variations in geomagnetic activity found in Jacchia '70 was
%erroneous, and is therefore not included here.
Te=Tg+Tl;

%Calculate temperature at desired altitude
Tx=371.668 + 0.0518806*Te - 294.3505*exp(-.00216222*Te);
A=2*(Te-Tx)/pi;
B=4.5e-6;
T1 = 1.9*(Tx-T0)/(zx-z0);
%T2 = 0 and so is not included.
T3 = -1.7*(Tx-T0)/((zx-z0)^3);
T4 = -.8*(Tx-T0)/((zx-z0)^4);
if Z>zx
    T=Tx + A*atan(T1*(Z-zx)*(1+B*(Z-zx)^2.5)/A);
else
    T=Tx + T1*(Z-zx) + T3*(Z-zx)^3 + T4*(Z-zx)^4;
end

%%% Calculate Standard Density
%Use bi-polynomial fit from Satellite Orbits
cij=getCoeff(Z,Te);
logRho=0;
for ii=0:5
    for jj=0:4
        logRho=logRho + cij(ii+1,jj+1)*((Z/1000)^ii)*((Te/1000)^jj);
    end
    
end

%Correction due to semi-annual density variation in thermosphere
Phi=(Jul1+Jul2 - 2400000.5 -36204)/365.2422; %number of tropical years since Jan 1, 1958
tauSA=Phi + 0.09544*((0.5+0.5*sin(2*pi*Phi+6.035))^1.65 - 0.5);

fZ=(5.876e-7*Z.^2.331 + 0.06328).*exp(Z*-2.868e-3);
gt=0.02835 + 0.3817*(1+0.4671*sin(2*pi*tauSA+4.137))*sin(4*pi*tauSA+4.259);

logRhoSA=fZ*gt;

%Correction due to geomagnetic activities
if Z<350
    logRhoGM=(0.012*Kp + 1.2e-5*exp(Kp))*(1-f);%Transition function from temp calc
else
    logRhoGM=0;
end

%Seasonal-latitude correction
logRhoSL=0.014*(Z-90).*exp(-0.0013*(Z-90).^2)*sin(1*pi*Phi - 1.72)*sin(lat)^3/abs(sin(lat));
    
logRho=logRho+logRhoSA+logRhoGM+logRhoSL;
rho=10^logRho;

%%% Pressure
% calculate mean molecular mass
M=cn(1);
for n=2:7
    M=M+cn(n)*(Z-90)^(n-1);
end

P=rho*Constants.molarGasConstant*T/M;
end

function c=getCoeff(Z,T)
%Retrieve coefficients for bi-polynomial fit taken from Tables 3.9 and 3.10
%of Satellite Orbits. These numbers were entered manually and may include
%some typing errors. This should be replaced by the process described in
%E. Gill, "Smooth Bi-Polynomial Interpolation of Jacchia 1971 Atmospheric
%Densities For Efficient Satellite Drag Computation," DLR-GSOC IB 96-1,
%German Aerospace Center (DLR), 1996
%This paper is unavailable at this time.

if Z<90 || Z>2500
    error('Z must be in range 90km < Z < 2500km')
end
if T<500 || T>1900
    error('T must be in range 500K < T < 1900K')
end

if T<850
    if Z<1000
        if Z<500
            if Z<180
                c=[-0.3520856e2  0.3912622e1 -0.8649259e2  0.1504119e3 -0.7109428e2
                      0.1129210e4  0.1198158e4  0.8633794e3 -0.3577091e4  0.1970558e4
                     -0.1527475e5 -0.3558481e5  0.1899243e5  0.2508241e5 -0.1968253e5
                      0.9302042e5  0.3646554e6 -0.3290364e6 -0.1209631e5  0.8438137e5
                     -0.2734394e6 -0.1576097e7  0.1685831e7 -0.4282943e6 -0.1345593e6
                      0.3149696e6  0.2487723e7 -0.2899124e7  0.1111904e7  0.3294095e4];
            else
                c=[ 0.2311910e2  0.1355298e3 -0.8424310e3  0.1287331e4 -0.6181209e3
                     -0.1057776e4  0.6087973e3  0.8690566e4 -0.1715922e5  0.9052671e4
                      0.1177230e5 -0.3164132e5 -0.1076323e4  0.6302629e5 -0.4312459e5
                     -0.5827663e5  0.2188167e6 -0.2422912e6  0.2461286e5  0.6044096e5
                      0.1254589e6 -0.5434710e6  0.8123016e6 -0.4490438e6  0.5007458e5
                     -0.9452922e5  0.4408026e6 -0.7379410e6  0.5095273e6 -0.1154192e6];
            end
        else
            c=[-0.1815722e4  0.9792972e4 -0.1831374e5  0.1385255e5 -0.3451234e4
                  0.9851221e4 -0.5397525e5  0.9993169e5 -0.7259456e5  0.1622553e5
                 -0.1822932e5  0.1002430e6 -0.1784481e6  0.1145178e6 -0.1641934e5
                  0.1298113e5 -0.7113430e5  0.1106375e6 -0.3825777e5 -0.1666915e5
                 -0.1533510e4  0.7815537e4  0.7037562e4 -0.4674636e5  0.3516946e5
                 -0.1263680e4  0.7265792e4 -0.2092909e5  0.2936094e5 -0.1491676e5];
        end
    else
        c=[ 0.3548698e3 -0.2508685e4  0.6252742e4 -0.6755376e4  0.2675763e4
             -0.5370852e3  0.4182586e4 -0.1151114e5  0.1338915e5 -0.5610580e4
             -0.2349586e2 -0.8941841e3  0.4417927e4 -0.6732817e4  0.3312608e4
              0.3407073e3 -0.1531588e4  0.2179045e4 -0.8841341e3 -0.1369769e3
             -0.1698470e3  0.8985697e3 -0.1704797e4  0.1363098e4 -0.3812417e3
              0.2494943e2 -0.1389618e3  0.2820058e3 -0.2472862e3  0.7896439e2];
    end
else
    if Z<1000
        if Z<500
            if Z<180
                c=[-0.5335412e2  0.2900557e2 -0.2046439e2  0.7977149e1 -0.1335853e1
                      0.1977533e4 -0.7091478e3  0.4398538e3 -0.1568720e3  0.2615466e2
                     -0.2993620e5  0.5187286e4 -0.1989795e4  0.3643166e3 -0.5700669e2
                      0.2112068e6 -0.4483029e4 -0.1349971e5  0.9510012e4 -0.1653725e4
                     -0.7209722e6 -0.7684101e5  0.1256236e6 -0.6805699e5  0.1181257e5
                      0.9625966e6  0.2123127e6 -0.2622793e6  0.1337130e6 -0.2329995e5];
            else
                c=[ 0.4041761e2 -0.1305719e3  0.1466809e3 -0.7120296e2  0.1269605e2
                     -0.8127720e3  0.2273565e4 -0.2577261e4  0.1259045e4 -0.2254978e3
                      0.5130043e4 -0.1501308e5  0.1717142e5 -0.8441698e4  0.1518796e4
                     -0.1600170e5  0.4770469e5 -0.5473492e5  0.2699668e5 -0.4870306e4
                      0.2384718e5 -0.7199064e5  0.8284653e5 -0.4098358e5  0.7411926e4
                     -0.1363104e5  0.4153499e5 -0.4793581e5  0.2377854e5 -0.4310233e4];
            end
        else
            c=[-0.4021335e2 -0.1326983e3  0.3778864e3 -0.2808660e3  0.6513531e2
                  0.4255789e3  0.3528126e3 -0.2077888e4  0.1726543e4 -0.4191477e3
                 -0.1821662e4  0.7905357e3  0.3934271e4 -0.3969334e4  0.1027991e4
                  0.3070231e4 -0.2941540e4 -0.3276639e4  0.4420217e4 -0.1230778e4
                 -0.2196848e4  0.2585118e4  0.1382776e4 -0.2533006e4  0.7451387e3
                  0.5494959e3 -0.6604225e3 -0.3328077e3  0.6335703e3 -0.1879812e3];
        end
    else
        c=[ 0.1281061e2 -0.3389179e3  0.6861935e3 -0.4667627e3  0.1029662e3
              0.2024251e3  0.1668302e3 -0.1147876e4  0.9918940e3 -0.2430215e3
             -0.5750743e3  0.8259823e3  0.2329832e3 -0.6503359e3  0.1997989e3
              0.5106207e3 -0.1032012e4  0.4851874e3  0.8214097e2 -0.6527048e2
             -0.1898953e3  0.4347501e3 -0.2986011e3  0.5423180e2  0.5039459e1
              0.2569577e2 -0.6282710e2  0.4971077e2 -0.1404385e2  0.8450500e0];
    end
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
