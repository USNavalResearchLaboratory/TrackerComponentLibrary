function e=seawaterRelPermittivity(T,v,S)
%%SEAWATERRELPERMITTIVITY Compute the relative permittivity (which is the
%           same as the permittivity divided by the permittivity of a
%           vacuum) of sea water as a function of temperature, the
%           electromagnetic frequency considered, and salinity. For nonzero
%           salinities, the approximation used here is valid from -2
%           degrees centigrade to 29 degrees centigrade, for frequencies up
%           to 90GHz. The model was fit to salinities from 20ppt to 40ppt,
%           but also matches the corresponding freshwater model given in
%           [1] when S=0. The sign convention used here is that the imag(e)
%           is negative, corresponding to an absorptive loss in the medium
%           (positive would be an amplification).
%
%INPUTS: T The temperature of the water in degrees Centigrade. The default
%          if omitted or an empty matrix is passed is 15.
%        v The frequency of the light given in Hertz. The default if
%          omitted or an empty matrix is passed is 2e9 (2 GHz).
%        S The salinity of the water given in parts per thousand (by
%          weight)of NaCl. The default if omitted or an empty matrix is
%          passed is 35. That corresponds to 35 grams of salt per Liter of
%          seawater.
%
%OUTPUTS: e The estimated relative permittivity. This is a unitless
%           quantity.
%
%The formulae in Section IV of [1] are implemented. The salinity fit is
%valid for S=0 (it simplifies to the freshwater solution in the paper) and
%for seawater, it was fit to salinities from 20ppt to 40ppt.
%
%The square root of the product of the relative permittivity and the
%relative permeability is the index of refraction. Due to the accuracy
%limitations in the models here and the fact that the relative pereability
%of water is close to 1, one can approximate the refractive index of water
%as just the square root of the permittivity.
%
%EXAMPLE:
%This plots the real component of the permittivity at 15 degrees  with
%20ppm salinity over a wide range of frequencies.
% numPts=1000;
% f=logspace(1,12,numPts);
% T=15;
% S=30;
% n=zeros(numPts,1);
% for k=1:numPts
%     n(k)=seawaterRelPermittivity(T,f(k),S);
% end
% figure(1)
% clf
% loglog(f,real(n),'-b','linewidth',4)
% xlabel('Frequency')
% ylabel('real(n)')
%
%REFERENCES:
%[1] T. Meissner and F. J. Wentz,"The complex dielectric constant of pure
%    water from microwave satellite observations," IEEE Transactions on
%    Geoscience and Remote Sensing, vol. 42, no. 9, pp. 1836-1849, Sep.
%    2004.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1||isempty(T))
   T=15;%Centigrade
end

if(nargin<2||isempty(v))
    v=2e9;%GHz
end

if(nargin<3||isempty(S))
    S=35;%ppt
end

%Convert from Hz to GHz.
v=v/1e9;

%Equation 12
sigma35=2.903602+8.607e-2*T+4.738817e-4*T^2-2.991e-6*T^3+4.3047e-9*T^4;

%Equation 15
alpha0=(6.9431+3.2841*S-9.9486e-2*S^2)/(84.850+69.024*S+S^2);
%Equation 16
alpha1=49.843-0.2276*S+0.198e-2*S^2;

%Equation 14
RT15Ratio=1+alpha0*(T-15)/(alpha1+T);

%Equation 13
R15=S*(37.5109+5.45216*S+1.4409e-2*S^2)/(1004.75+182.283*S+S^2);

%Equation 11.
sigma=sigma35*R15*RT15Ratio;

%Table III in [1].
a0=5.7230;
a1=2.2379e-2;
a2=-7.1237e-4;
a3=5.0478;
a4=-7.0315e-2;
a5=6.0059e-4;
a6=3.6143;
a7=2.8841e-2;
a8=1.3652e-1;
a9=1.4825e-3;
a10=2.4166e-4;

%Equation 7 in [1] for freshwater.
eS=(3.70886e4-8.2168e1*T)/(4.21854e2+T);
%Equation 8 in [1] for freshwayer.
e1=a0+a1*T+a2*T^2;
v1=(45+T)/(a3+a4*T+a5*T^2);
eInf=a6+a7*T;
v2=(45+T)/(a8+a9*T+a10*T^2);

%Table VI in [1].
b0=-3.56417e-3;
b1=4.74868e-6;
b2=1.15574e-5;
b3=2.39357e-3;
b4=-3.13530e-5;
b5=2.52477e-7;
b6=-6.28908e-3;
b7=1.76032e-4;
b8=-9.22144e-5;
b9=-1.99723e-2;
b10=1.81176e-4;
b11=-2.04265e-3;
b12=1.57883e-4;

%Equation 17 in [1], scaling the freshwater terms for saltwater.
eS=eS*exp(b0*S+b1*S^2+b2*T*S);
v1=v1*(1+S*(b3+b4*T+b5*T^2));
e1=e1*exp(b6*S+b7*S^2+b8*T*S);
v2=v2*(1+S*(b9+b10*T));
eInf=eInf*(1+S*(b11+b12*T));

%The value of 1/(2*pi*e0) used in [1] (see the paragraph after
%Equation 3). e0 is the permittivity of a vacuum.
invTwoPie0=17.97510;%(GHz m/S).

%Equation 6 in 1. The sigma term is 0 for freshwater.
e=(eS-e1)/(1+1i*v/v1)+(e1-eInf)/(1+1i*v/v2)+eInf-1i*sigma/(invTwoPie0*v);

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
