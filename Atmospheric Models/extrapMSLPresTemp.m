function [P0,T0]=extrapMSLPresTemp(P,z,T)
%%EXTRAPMSLPRESTEMP  Given a measured pressure P at orthometric height z,
%         extrapolate the pressure at mean sea level (MSL) P0, which is
%         orthometric height 0. If the temperature at height z is known,
%         then the temperature T0 at MSL can also be extrapolated. This
%         uses the U.S. Standard Atmosphere 1976 and is valid up to 86km.
%         If a temperature at altitude is not given, then it is assumed
%         that T0 is the standard temperature in Constants.standardTemp.
%
%INPUTS: P A pressure measurement in Pascals at height z.
%        z An orthometric height in meters where the pressure measurement
%          was taken.
%        T An optional temperature measurement in Kelvin at height z. If
%          omitted, it is assumed that the temperature at MSL is
%          Constants.standardTemp.
%
%OUTPUTS: P0 The air pressure at MSL extrapolated using the U.S. Standard
%            Atmosphere 1976.
%         T0 The extrapolated temperature at MSL. If T is not provided,
%            then T0 will just be Constants.standardTemp
%
%%The atmospheric model is defined in [1], and the full model goes up to
%1000km altitude. Up to an altitude of 32km, the U.S. standard atmosphere
%is identical to the International Civil Aviation Organization's (ICAO)
%standard atmosphere.
%
%REFERENCES:
%[1] U.S. Standard Atmosphere, 1976, National Oceanic and Atmospheric
%    Administration Std., Oct. 1976.
%
%June 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    T=Constants.standardTemp;
end

if(z>86e3)
   warning('The model has not been implemented for heights above 86km') 
end

%%%Constants defined in the standard atmosphere.
%The radius of the spherical Earth used in the standard atmosphere. THis is
%not the same as the semi-major axis of the WGS-84 reference ellipsoid.
rE=6356766;%Meters
%The 1962 estimate of the UNiversal Gas Constant that is used in the
%standard.
R=8.31432;%J/(mole*K)
%The mean moleacular mass of air at sea level.
M0=28.9644e-3;%kg/mole
%Standard acceleration due to gravity.
g0=9.80665;%m/s

%Hb and LMb are from table 4 in the U.S. Standard Atmosphere 1976.
%Hb is the base geopotential height and LMb is the molecular scale
%temperature gradient from that base.
HbTable=[0;
    11000;
    20000;
    32000;
    47000;
    51000;
    71000;
    84852];
LMbTable=[-6.5e-3;
    0;
    1e-3;
    2.8e-3;
    0;
    -2.8e-3;
    -2e-3;
    0];

%Get the geopotential height using the approximation in the standard.
H=z*(rE/(rE+z));

%Next, we have to determine the height Hb of the base level.
HIdx=7;
for n=2:7
    if(H<HbTable(n))
       HIdx=n-1;
       break;
    end
end

%Next, fill in a whole TMb table up to geopotential height Hb. If T is
%given, then we must go backwards to T0. If T is not given, then we must
%start with T0 and go forwards. This uses equation 23 in the standard.
TMbTable=zeros(HIdx,1);
if(nargin>2)%If T is given.
    TMbTable(HIdx)=T-LMbTable(HIdx)*(H-HbTable(HIdx));
    for n=(HIdx-1):-1:1
        TMbTable(n)=TMbTable(n+1)-LMbTable(n)*(HbTable(n+1)-HbTable(n));
    end
else%If T is not given, start with the standard temperature and go forward.
    TMbTable(1)=Constants.standardTemp;
    for n=2:8
        TMbTable(n)=TMbTable(n-1)+LMbTable(n-1)*(HbTable(n)-HbTable(n-1));
    end
end

T0=TMbTable(1);

%Go from the pressure P at H to the pressure Pb at Hb and then loop
%backwards to P0.
if(LMbTable(HIdx)==0)
    P=P/exp(-g0*M0*(H-HbTable(HIdx))/(R*TMbTable(HIdx)));
else
    P=P/(TMbTable(HIdx)/(TMbTable(HIdx)+LMbTable(HIdx)*(H-HbTable(HIdx))))^(g0*M0/(R*LMbTable(HIdx)));
end

for n=(HIdx-1):-1:1
    if(LMbTable(n)==0)
        P=P/exp(-g0*M0*(HbTable(n+1)-HbTable(n))/(R*TMbTable(n)));
    else
        P=P/(TMbTable(n)/(TMbTable(n)+LMbTable(n)*(HbTable(n+1)-HbTable(n))))^(g0*M0/(R*LMbTable(n)));
    end
end

P0=P;
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
