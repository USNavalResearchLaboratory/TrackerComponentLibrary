function [N,rho,P]=atmosParam4GasTemp(gasTable,T,rhow)
%%ATMOSPARAM4GASTEMP  Get basic parameters for atmospheric refractivity,
%                     density, and pressure given a table of constitudent
%                     gasses in the atmosphere, the temperature and the
%                     absolute humidity. The model is best suited for
%                     altitudes below 90km as ionozed parameters, such as
%                     anomolous oxygen, are not used.
%
%INPUTS: gasTable An NX2 cell array where gasTable{i,1} is a string
%                 describing the ith constituent atmospheric element and
%                 gasTable{i,2} is the number density of the element in
%                 particles per cubic meter. For a list of constituent
%                 elements that can be handled, see the documentation for
%                 the Constants.gasProp method. Unknown constituent
%                 elements that are passed will be ignored.
%               T The temperature in degrees Kelvin.
%            rhow An optional parameter specifying the mass density of
%                 water vapor at the point in question in units of
%                 kilograms per cubic meter. If omitted, the air is assumed
%                 to be dry (rhow=0). The total density of the air is
%                 assumed to be the sum of the dry air density and rhow.
%                 Alternatively, this parameter can be omitted and 'H2O'
%                 can be given as one of the constituent elements in
%                 gasTable.
%
%OUTPUTS: N The refractivity of the atmosphere. In this model, N is always
%           real. N=10^6*(n-1) where n is the index of refraction. This is
%           generally valid for frequencies from L-band (1GHz) to 10 GHz
%           (the middle of X-band).
%       rho The atmospheric density at the point in question in units of
%           kilograms per cubic meter. 
%         P The atmospheric pressure at the point in question in units of
%           Newtons per square meter (Pascals). It assumes that the gasses
%           can be treated as ideal gasses.
%
%The refractive index is then found using the dry and wet air densities
%using the formula of [1], which should be valid for frequencies between
%1GHz and 10GHz. It ignores the lossiness of the atmosphere.
%
%REFERENCES:
%[1] J. M. Aparicio and S. Laroche, "An evaluation of the expression of
%    the atmospheric refractivity for GPS signals," Journal of Geophysical
%    Research, vol. 116, no. D11, 16 Jun. 2011.
%
%September 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
   rhow=0; 
end

numGasses=size(gasTable,1);
numberDensities=[gasTable{:,2}];
totalGasNumberDensity=sum(numberDensities);

%rhod will be the total mass density of dry air at the point in kg/m^3.
rhod=0;
for curGas=1:numGasses
    AMU=Constants.gasProp(gasTable{curGas,1});
    %If an unknown gas is given, then ignore it.
    if(~isempty(AMU))
        %The number density has units of particles per cubic meter.
        %Multiplied by the atomic mass of the gas, we have atomic mass
        %units (AMU) per cubic meter. Multiplied by the value of the atomic
        %mass unit in kilograms, we get kilograms per cubic meter.
        rhod=rhod+Constants.atomicMassUnit*numberDensities(curGas)*AMU;
    else
        %This is so that the number density of the unknown constitutent is
        %also ignored when computing the pressure, so that the results are
        %consistent with the computation of the density.
        numberDensities(curGas)=0;
    end
end

rho=rhod+rhow;

%To use the ideal gas law to find the air pressure, the number of water
%molecules per cubic meter of the atmosphere is needed. This is obtained
%using the molar mass of water (H2O) and Avogadro's constant
Av=Constants.AvogadroConstant;
%The molar mass of water (H2O) in atomic mass units (grams per mole).
HAMU=Constants.elementAMU(1);
OAMU=Constants.elementAMU(8);
MMWater=HAMU*2+OAMU;

%The number of atoms of water per cubic meter. The 1000 transforms grams to
%kilograms.
NH2O=rhow/(1000*Av*MMWater);

%The total number density of the gasses in the atmosphere. That is, the
%number of atoms per cubic meter.
NTotal=totalGasNumberDensity+NH2O;
kB=Constants.BoltzmannConstant;
P=NTotal*kB*T;

%T is the temperature in Kelvin.
tau=273.15/T-1;
N0=(222.682+0.069*tau)*rhod+(6701.605+6385.886*tau)*rhow;
N=N0*(1+10^(-6)*N0/6);
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
