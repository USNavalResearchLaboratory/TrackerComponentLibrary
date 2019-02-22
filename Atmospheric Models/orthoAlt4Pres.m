function z=orthoAlt4Pres(P,P0,T0)
%%ORTHOALT4PRES Get the orthometric height (pressure altitude) for a
%               point up to 86km altitude using the U.S. Standard
%               Atmosphere 1976 model given a reading of the air pressure
%               at the point. The model is built into many altimeters.
%
%INPUTS: P A vector or matrix of the pressures in Pascals that are to be
%          converted into orthometric heights in meters.
%       P0 An optional parameter specifying the pressure at sea level at
%          the location in question in Pascals. This pressure might have 
%          to be interpolated down below terrain given a measurement on
%          the terrain. If omitted, the value of
%          Constants.standardAtmosphericPressure is used for the pressure.
%          Except in certain instances, aircraft over 18,000ft altitude in
%          the United States have to have their altimeter calibrated to
%          the standard atmospheric pressure and not the true pressure on
%          the ground/ at sea level.
%       T0 An optional parameter specifying the temperature at sea level
%          in degrees Kelvin. The standard atmosphere only specifies one
%          value, which is in Constants.standardTemp and is the default
%          value used if this parameter is omitted.
%
%OUTPUTS: z The orthometric height of the point implied by the specified air
%           pressure.The standard does not specify a reference geoid/
%           permanent tide system.
%
%The atmospheric model is defined in [1]. and the full model goes up to
%1000km altitude. Up to an altitude of 32km, the U.S. standard atmosphere
%is identical to the International Civil Aviation Organization's (ICAO)
%standard atmosphere. Regulations related to how altimeters should be set
%with respect to a standard pressure are at
%https://www.faa.gov/air_traffic/publications/atpubs/aim/aim0702.html
%
%REFERENCES:
%[1] U.S. Standard Atmosphere, 1976, National Oceanic and Atmospheric
%   Administration Std., Oct. 1976.
%
%June 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2)
    P0=Constants.standardAtmosphericPressure;%N/m
end

if(nargin<3)
    T0=Constants.standardTemp;
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

numP=size(P(:));
z=zeros(size(P));
for curP=1:numP
    PCur=P(curP);

    %To determine the correct formula to use, we must first determine which
    %band of altitudes the pressure falls. Thus, we will first evaluate
    %pres4OrthoAlt for the pressure at all of the bands that separate, knowing
    %that the pressure monotonically decreases with altitude in the model.
    H=11000;
    zt=rE*(H/(rE-H));
    [Ph,~,Hb,Pb,TMb,LMb]=presTemp4OrthoAlt(zt,P0,T0);
    if(PCur<Ph)
        H=20000;
        zt=rE*(H/(rE-H));
        [Ph,~,Hb,Pb,TMb,LMb]=presTemp4OrthoAlt(zt,P0,T0);
        if(PCur<Ph)
            H=32000;
            zt=rE*(H/(rE-H));
            [Ph,~,Hb,Pb,TMb,LMb]=presTemp4OrthoAlt(zt,P0,T0);
            if(PCur<Ph)
                H=47000;
                zt=rE*(H/(rE-H));
                [Ph,~,Hb,Pb,TMb,LMb]=presTemp4OrthoAlt(zt,P0,T0);
                if(PCur<Ph)
                    H=51000;
                    zt=rE*(H/(rE-H));
                    [Ph,~,Hb,Pb,TMb,LMb]=presTemp4OrthoAlt(zt,P0,T0);
                    if(PCur<Ph)
                        H=71000;
                        zt=rE*(H/(rE-H));
                        [Ph,~,Hb,Pb,TMb,LMb]=presTemp4OrthoAlt(zt,P0,T0);
                        if(PCur<Ph)
                            H=72000;
                            zt=rE*(H/(rE-H));
                            [~,~,Hb,Pb,TMb,LMb]=presTemp4OrthoAlt(zt,P0,T0);
                        end
                    end
                end
            end
        end
    end

    %Solve the inverse of the appropriate air pressure equation to determine
    %the geopotential height.
    if(LMb==0)
        %The inverse of the exponential expression for air pressure as a
        %function of altitude in Equation 33b
        H=Hb-R*TMb*log(PCur/Pb)/(g0*M0);
    else
        %The inverse of air pressure according to equation 33a, which is the
        %linear formulation. 
        H=Hb+(TMb/LMb)*((PCur/Pb)^(-LMb*R/(g0*M0))-1);
    end

    %Convert the geopotential height into an orthometric height using the
    %approximation in the standard.
    z(curP)=rE*(H/(rE-H));
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
