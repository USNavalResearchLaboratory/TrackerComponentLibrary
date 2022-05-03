function c=speedOfSoundInAir(algorithm,T,input3,input4,CO2Frac)
%%SPEEDOFSOUND Compute the speed of sound for given atmospheric parameters
%              using the selected algorithm.
%
%INPUTS:    If no inputs are given, then a value for the speed of sound at
%           standard temperature and pressure from [1] is returned.
%           Otherwise, the inputs are:
%           algorithm This is a number from 0 to 2 specifying which
%                     algorithm is to be used. References for the
%                     algorithms are given below. The algorithms are
%                     0  A detailed approach utilizing the temperature,
%                        pressure and the gas content of the air.
%                     1  A simple ideal-gas approximation utilizing the
%                        temperature and relative humidity and assuming
%                        standard pressure (101325 Pascals).
%                     2  An approximation to the method of algorithm 0.
%                     All of the algorithms have limitations on their range
%                     of validity, as described below.
%                   T The temperature in degrees Kelvin.
%              input3 This input depends on the algorithm chosen. For
%                     each algorithm this is respectively:
%                     0 P The pressure in Pascals (kg/(m*s^2));
%                     1 relHumid The relative humidity as a value from
%                        0-1.
%                     2 P The pressure in Pascals (kg/(m*s^2)).
%              input4 This input depends on the algorithm chosen. For
%                     each algorithm this is respectively:
%                     0 gasTable An NX2 cell array where gasTable{i,1} is a
%                       string describing the ith constituent atmospheric
%                       element and gasTable{i,2} is the number density of
%                       the element in particles per cubic meter. Note that
%                       only the relative number densities matter, so
%                       multiplying all of the number densities by a
%                       constant does not change the final result. For a
%                       list of constituent elements that can be handled,
%                       see the documentation for the Constants.gasProp
%                       method. See the sample code in the comments below
%                       for an example of how this input is used. Unknown
%                       constituent elements that are passed will be
%                       ignored.
%                     1 This input is not used.
%                     2 H2OFrac The water vapor mole fraction. This is the
%                       fraction of H2O molecules per cubic meter to total
%                       molecules per cubic meter of air. If this parameter
%                       is omitted, then a dry (H2OFrac=0) atmosphere is
%                       used.
%             CO2Frac This input is only used by algorithm 2. This is 
%                     the carbon dioxide mole fraction, which is the
%                     fraction of CO2 molecules per cubic meter to total
%                     molecules per cubic meter of air. If this parameter
%                     is omitted, then a CO2-free atmosphere (CO2Frac=0) is
%                     used.
%
%OUTPUTS: c The speed of sound under the specified conditions in meters per
%           second.
%
%If no inputs are provided, then the value for the speed of sound is the
%value taken at standard temperature and pressure (0 degrees Centigrade,
%no humidity, 1 atmosphere [101325 Pascals]) as in [1].
%
%If algorithm 0 is chosen, then the solution based upon individual gas
%concentrations at 1 atmosphere pressure from [2] is used where the
%necessary physical constants are taken from the Constants class.
%
%If algorithm 1 is chosen, the ideal gas approximation, then the method of
%[3] is used. The algorithm assumes a standard atmospheric pressure of
%101325Pa.
%
%If algorithm 2 is chosen, then the simple polynomial approximation of the
%[2] is used.
%
%Algorithm 0 is valid for most combinations of parameters, but is limited
%by the tabulated values of the second virial coefficient and the
%ideal gas specific heat in the Constants class for various gasses. Those
%values are often not available for certain constituent gasses at
%temperatures deviating greatly from the freezing point (273.15K),
%particularly below freezing.
%Algorithm 1 was derived for temperatures from 273.15K (0C) to 303.15K
%(30C) and assumes standard pressure (101325 Pa).
%Algorithm 2 was derived for temperatures from 273.15K (0C) to 303.15K
%(30C), pressures from 75000Pa to 102000Pa, water vapor mole fractions up
%to 0.06 (6%) and carbon dioxide mole fractions up to 0.01 (1%).
%
%The method can be used as
%c=speedOfSoundInAir();
%or
%c=speedOfSoundInAir(0,T,P,gasTable);
%or
%c=speedOfSoundInAir(1,T,relHumid);
%or
%c=speedOfSoundInAir(2,T,P,H2OFrac,CO2Frac);
%
%EXAMPLE: Suppose that we would like to provide a table of constituent
%gasses based on a standard atmospheric model (The NRLMSISE-00 model).
%Given a temperature, the number densities implied by the model imply a
%pressure. However, we would like to use a measured local pressure. Also,
%we have a local measurement of the humidity (the NRLMSISE-00 model assumes
%zero humidity). The model depends on the observer's location and the time.
%For example:
% %0:00 UTC on 1 January 2000 (Gregorian calendar)
% year=2000;
% month=1;
% day=1;
% hour=0;
% minute=0;
% second=0;
% [Jul1,Jul2]=Cal2UTC(year,month,day,hour,minute,second);
% [~,dayCount,secondInDay]=UTC2DayCount(Jul1,Jul2);
% 
% %The latitude and longitude of the Mauna Kea Observatory in Hilo, Hawaii.
% phi=19.823*pi/180;%North latitude in radians.
% lambda=-155.470*pi/180;%East longitude in radians.
% height=4200;%Assume around 4.2km ellipsoidal height
% 
% %Get the standard atmospheric parameters for the given location and time.
% [gasTable,t,d]=NRLMSISE00GasTemp(dayCount,secondInDay,[phi;lambda;height]);
% 
% %Let's assume that the temperature at the observer is 10 degrees C,
% %rather than the standard temperature that could have been returned by
% %the NRLMSISE00GasTemp function.
% TC=10;%10 degrees centigrade
% T=convertTemperatureUnits(TC,'C','K');%Convert to Kelvin
% 
% %The standard atmospheric parameters assume zero humidity. Let's assume
% %that we measured the humidity to be 10%. We need the number density of
% %the H2O in the air.
% relHumid=0.1;
% absHumid=relHumid2AbsHumid(relHumid,T);
% numberDensH2O=absHumid2NumberDensH2O(absHumid);
% %We will add the water to the gas table.
% gasTable{end+1,1}='H2O';
% gasTable{end,2}=numberDensH2O;
% 
% %We need the air pressure. The pressure implied by the temperatue and the
% %table of gasses is
% [~,~,P]=atmosParam4GasTemp(gasTable,T);
% %However, note that since the speedOfSoundInAir function only depends on
% %the relative abundance of the individual gasses, one could use a
% %different, measured air pressure, assuming that the relative abundance
% %has not changed.
% %Finally, we can now find the speed of sound in air
% c=speedOfSoundInAir(0,T,P,gasTable)
%
%REFERENCES:
%[1] D. H. Smith and R. G. Harlow, "The velocity of sound in air, nitrogen
%    and argon," British Journal of Applied Physics, vol. 14, no. 2, pp.
%    102-106, Feb. 1963.
%[2] O. Cramer, "The variation of the specific heat ratio and the speed of
%    sound in air with temperature, pressure, humidity, and CO2
%    concentration," Journal of the Acoustical Society of America, vol. 93,
%    no. 5, pp. 2510-2516, May 1993.
%[3] G. S. K. Wong and T. F. W. Embleton, "Variation of the speed of sound
%    in air with humidity and temperature," Journal of the Acoustical
%    Society of America, vol. 77, no. 5, pp. 1710-1712, May 1985.
%
%March 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Use the value for STP from Smith and Harlow if nothing else is provided.
if(nargin<1)
    c=331.45;
    return;
end

switch(algorithm)
    case 0
        P=input3;
        gasTable=input4;
        numGasses=size(gasTable,1);
        
        A=Constants.AvogadroConstant;

        %Compute the mass density and total mass density in kg/m^3
        massDensityList=zeros(numGasses,1);
        for curGas=1:numGasses
            %The atomic weight in amu (kg/mol) of the constituent gas. If
            %the gas is not in the Constants class, then it is not included
            %in the speed of sound computation.
            atomWeight=Constants.gasProp(gasTable{curGas,1})/1000;
            
            if(~isempty(atomWeight))
                %Convert particles/m^3 to moles to m^3 and then multiply by
                %the atomic weight to get the mass density in kg/m^3.
                massDensityList(curGas)=(gasTable{curGas,2}/A)*atomWeight;
            end
        end
        totalMassDensity=sum(massDensityList);
        
        c2=0;
        for curGas=1:numGasses
            gasName=gasTable{curGas,1};
            massDensity=massDensityList(curGas);
            [M,C0p,B,dBdT,d2BdT2]=Constants.gasProp(gasName,T);
            %M is in amu (g/mol). C0p is in J/(kg*K). B is in m^3/mol.
            
            %If values for the constituent element are tabulated.
            if(~isempty(M))
                %Convert the ideal gas specific heat from J/(kg*K) to
                %kJ/(kg*K);
                C0p=C0p/1e3;

                %The ideal gas constant in J/(mol*K)
                R=Constants.molarGasConstant;

                %Equation 5 (kJ/(kg*K))
                C1p=C0p-(R/M)*(P/(R*T))*T^2*d2BdT2;
                %Equation 6 (kJ/(kg*K))
                C1v=C1p-(R/M)*(1+2*P/(R*T)*T*dBdT);

                %The unitless specific heat ratio.
                gamma=C1p/C1v;

                %Convert from AMU (g/mol) to kg/mol.
                M=M/1e3;

                %Equation 8, summed weighted for each gas based on its relative
                %mass density
                c2=c2+(massDensity/totalMassDensity)*gamma*(R*T/M)*(1+2*P*B/(R*T));
            end
        end
        
        c=sqrt(c2);
    case 1%Use the ideal gas approximation of the Wong and Embleton paper.
        relHumid=input3;
        
        %Convert T from degrees Kelvin to degrees Centigrade
        t=T+Constants.absoluteZero;
        
        if(t<0||t>30)
           warning('The temperature supplied is outside of the range used (0-30 degrees C) in the paper deriving the ideal gas approximation. The results might have reduced accuracy.')
        end
        
        %Equation 4
        At=9.2e-5+5.5e-6*t+4.25e-7*t^2;
        
        %Equation 3, units of moles per gram
        gammahdMh=0.04833+(relHumid-0.023*At);
        
        %Convert to moles per kilogram.
        gammahdMh=gammahdMh*1000;
        
        %The ideal gas constant in J/(mol*K)
        R=Constants.molarGasConstant;
        
        %Equation 1
        c=sqrt(gammahdMh*R*T);
        return;
    case 2%Use the approximate algorithm of the Cramer paper.
        P=input3;
        H2OFrac=input4;
        
        %Convert T from degrees Kelvin to degrees Centigrade
        t=T+Constants.absoluteZero;
        if(t<0||t>30)
           warning('The temperature supplied is outside of the range used (0-30 degrees C) in the paper deriving speed approximation. The results might have reduced accuracy.')
        end

        if(P<75000||P>102000)
            warning('The pressure supplied is outside of the range used (75000-102000 P) in the paper deriving the speed approximation. The results might have reduced accuracy.')
        end
        
        if(H2OFrac<0)
            error('Invalid mole fraction of water provided.');
        end
        
        if(H2OFrac>0.06)
            warning('The pressure supplied is outside of the range used (0-0.06) in the paper deriving the speed approximation. The results might have reduced accuracy.')
        end
        
        if(CO2Frac<0)
            error('Invalid mole fraction of carbon dioxide provided.');
        end
        
        if(CO2Frac>0.01)
            warning('The pressure supplied is outside of the range used (0-0.01) in the paper deriving the speed approximation. The results might have reduced accuracy.')
        end

    %From Table III for c
        a=[331.5024;
              0.603055;
             -0.000528;
             51.471935;
              0.1495874;
             -0.000782;
             -1.82e-7;
              3.73e-8;
             -2.93e-10;
             -85.20931;
             -0.228525;
              5.91e-5;
             -2.835149;
             -2.15e-13;
             29.179762;
              0.000486];

    %Equation 15
        c=a(1)+a(2)*t+a(3)*t^2+(a(4)+a(5)*t+a(6)*t^2)*H2OFrac...
            +(a(7)+a(8)*t+a(9)*t^2)*P+(a(10)+a(11)*t+a(12)*T^2)*CO2Frac...
            +a(13)*H2OFrac^2+a(14)*P^2+a(15)*CO2Frac^2+a(16)*H2OFrac*P*CO2Frac;
        return;
    otherwise
        error('Invalid algorithm number entered')
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
