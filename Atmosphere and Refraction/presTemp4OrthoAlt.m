function [P,T,Hb,Pb,TMb,LMb]=presTemp4OrthoAlt(z,P0,T0)
%%PRESTEMP4ORTHOALT Get the standard air pressure and temperature given an
%               orthometric altitude of a point up to 86km altitude
%               using the U.S. Standard Atmosphere 1976. The model is built
%               into many altimeters.
%
%INPUTS: z A vector or matrix of orthometric heights of the points (height
%          above mean sea level) in meters to convert to pressures. The
%          standard does not specify a reference geoid/ permanent tide
%          system.
%       P0 An optional parameter specifying the pressure at sea level at
%          the location in question in Pascals. This pressure might have 
%          to be interpolated down below terrain given a measurement on
%          the terrain. If omitted, the value of
%          Constants.standardAtmosphericPressure is used for the pressure.
%          Except in certain instances, aircraft over 18,000ft altitude in
%          the United States have to have their altimeter calibrated to
%          the standard atmospheric pressure and not the true pressure on
%          the ground/ at sea level.
%       T0 An optional parameter specifying the temperature at sea level in
%          degrees Kelvin. The standard atmosphere only specifies one
%          value, which is in Constants.standardTemp and is the default
%          value used if this parameter is omitted.
%
%OUTPUTS: P The standard pressure of the atmosphere at the specified
%           points, given in Pascals.
%         T The standard temperature in Kelvin at the specified altitude.
%        Hb The height of the base geopotential height used in the
%           computation. This is used in the function orthoAlt4Pres.
%        Pb The pressure at the geopotential height Hb. This is used in the
%           function orthoAlt4Pres
%       TMb The temperature used for the base geopotential height. This is
%           used in the function orthoAlt4Pres.
%       LMb The is the molecular scale temperature gradient from the base
%           used in the computation. This is used in the function
%           orthoAlt4Pres.
%
%The atmospheric model is defined in [1], and the full model goes up to
%1000km altitude. Up to an altitude of 32km, the U.S. standard atmosphere
%is identical to the International Civil Aviation Organization's (ICAO)
%standard atmosphere. Regulations related to how altimeters should be set
%with respect to a standard pressure are at
%https://www.faa.gov/air_traffic/publications/atpubs/aim/aim0702.html
%
%REFERENCES:
%[1] U.S. Standard Atmosphere, 1976, National Oceanic and Atmospheric
%    Administration Std., Oct. 1976.
%
%June 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    P0=Constants.standardAtmosphericPressure;%N/m
end

if(nargin<3)
    T0=Constants.standardTemp;
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

%TMb is computed using Equation 23 starting from T0. It isn't just
%tabulated in case the user wants to specify a non-standard T0.
TMbTable=zeros(8,1);
TMbTable(1)=T0;
for n=2:8
    TMbTable(n)=TMbTable(n-1)+LMbTable(n-1)*(HbTable(n)-HbTable(n-1));
end

%Compute the Pb values as describes on page 12
PbTable=zeros(8,1);
PbTable(1)=P0;
for n=2:8
    if(LMbTable(n-1)==0)
        PbTable(n)=PbTable(n-1)*exp(-g0*M0*(HbTable(n)-HbTable(n-1))/(R*TMbTable(n-1)));
    else
        PbTable(n)=PbTable(n-1)*(TMbTable(n-1)/(TMbTable(n-1)+LMbTable(n-1)*(HbTable(n)-HbTable(n-1))))^(g0*M0/(R*LMbTable(n-1)));
    end
end

numZ=length(z(:));
P=zeros(size(z));
T=zeros(size(z));
Hb=zeros(size(z));
Pb=zeros(size(z));
TMb=zeros(size(z));
LMb=zeros(size(z));
for curZ=1:numZ
    %Get the geopotential height using the approximation in the standard.
    H=z(curZ)*(rE/(rE+z(curZ)));
    
    %84852 is the geopotential height that is precisely the same as 86km
    %orthometric height, using the approximate conversion in the standard.
    %Note that the standard uses less than, not less than or equal to. However,
    %since there are no bifurcations, it does not matter. The reason for using
    %<= instead of < is that it makes it easy to find the correct pressure
    %region to use when implementing the inverse function orthoAlt4Pres.
    if(H<=11000)
        Pb(curZ)=PbTable(1);
        Hb(curZ)=HbTable(1);
        TMb(curZ)=TMbTable(1);
        LMb(curZ)=LMbTable(1);
    elseif(H<=20000)
        Pb(curZ)=PbTable(2);
        Hb(curZ)=HbTable(2);
        TMb(curZ)=TMbTable(2);
        LMb(curZ)=LMbTable(2);
    elseif(H<=32000)
        Pb(curZ)=PbTable(3);
        Hb(curZ)=HbTable(3);
        TMb(curZ)=TMbTable(3);
        LMb(curZ)=LMbTable(3);
    elseif(H<=47000)
        Pb(curZ)=PbTable(4);
        Hb(curZ)=HbTable(4);
        TMb(curZ)=TMbTable(4);
        LMb(curZ)=LMbTable(4);
    elseif(H<=51000)
        Pb(curZ)=PbTable(5);
        Hb(curZ)=HbTable(5);
        TMb(curZ)=TMbTable(5);
        LMb(curZ)=LMbTable(5);
    elseif(H<=71000)
        Pb(curZ)=PbTable(6);
        Hb(curZ)=HbTable(6);
        TMb(curZ)=TMbTable(6);
        LMb(curZ)=LMbTable(6);
    else
        Pb(curZ)=PbTable(7);
        Hb(curZ)=HbTable(7);
        TMb(curZ)=TMbTable(7);
        LMb(curZ)=LMbTable(7);
    end

    if(LMb(curZ)==0)
        %The exponential expression for air pressure as a function of
        %altitude. Equation 33b
        P(curZ)=Pb(curZ)*exp(-g0*M0*(H-Hb(curZ))/(R*TMb(curZ)));
    else
        %Air pressure according to equation 33a, which is the linear
        %formulation. 
        P(curZ)=Pb(curZ)*(TMb(curZ)/(TMb(curZ)+LMb(curZ)*(H-Hb(curZ))))^(g0*M0/(R*LMb(curZ)));
    end
    T(curZ)=TMb(curZ)+LMb(curZ)*(H-Hb(curZ));
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
