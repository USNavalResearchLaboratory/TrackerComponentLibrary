function e=freshwaterRelPermittivity(T,v,algorithm)
%%FRESHWATERPERMITTIVITY Compute the relative permittivity (which is the
%           same as the permittivity divided by the permittivity of a
%           vacuum) of freshwater as a function of temperature and the
%           electromagnetic frequency considered. The sign convention used
%           here is that the imag(e) is negative, corresponding to an
%           absorptive loss in the medium (positive would be an
%           amplification).
%
%INPUTS: T The temperature of the water in degrees Centigrade. The default
%          if omitted or an empty matrix is passed is 15. This input isn't
%          used if algorithm=2.
%        v The frequency of the light given in Hertz. The default if
%          omitted or an empty matrix is passed is 2e9 (2 GHz).
% algorithm This specifies which model to use. possible values are
%          0 (The default if omitted or an empty matrix is passed and
%            v<=25e12Hz) Use the approximation in [1], which is valid for
%            temperatures between -5 and 100 degrees centigrade and
%            frequencies from 0 to 25THz.
%          1 Use the approximation in [2], which is valid for temperatures
%            between -20 (supercooled water) and 40 degrees centigrade for
%            frequencies up to 500GHz. 
%          2 (The default if omitted or an empty matrix is passed and 
%            v>25e12 Hz) Use the approximation of [3], which does not rely
%            on temperature and goes from about 1e10Hz to 1e18Hz (above
%            ultraviolet).
%          3 Use the solution of [4]. Data, constants and modified
%            equations from the original authors (not in the paper) were
%            required to recreate the results. It actually computes sqrt(e)
%            (an approximation of the refractive index). When considering
%            sqrt(e), the complex part is a pretty good fit from 1e7 to
%           1e16Hz. The real part is an okay fit from 1e7 to to 1e14Hz and
%            then from about 1e16 to 1e18Hz. The model uses temperature,
%            but was not validated over any particular range. Compared to
%            other models, it appears to be  biased at low end of its
%            frequency range.
%
%OUTPUTS: e The estimated relative permittivity. This is a unitless
%           quantity.
%
%The square root of the product of the relative permittivity and the
%relative permeability is the index of refraction. Due to the accuracy
%limitations in the models here and the fact that the relative pereability
%of water is close to 1, one can approximate the refractive index of water
%as just the square root of the permittivity.
%
%EXAMPLE:
%This recreates the first plot in Figure 3 in [3] and also adds the line
%for algorithm 0 at 15 degrees in its valid range.
% numPts=1000;
% f=logspace(10,18,numPts);
% n0=zeros(numPts,1);
% n2=zeros(numPts,1);
% for k=1:numPts
%     n0(k)=freshwaterRelPermittivity(15,f(k),0);
%     n2(k)=freshwaterRelPermittivity([],f(k),2);
% end
% figure(1)
% clf
% semilogx(f,log10(real(n2)),'-r','linewidth',4)
% hold on
% %Only display the region where algorithm 0 is valid.
% sel=(f<=25e12);
% n0(~sel)=NaN;
% semilogx(f,log10(real(n0)),'-b','linewidth',2)
% xlabel('Frequency')
% ylabel('log_{10}(real(n))')
% legend('Algorithm 2', 'Algorithm 0')
%
%REFERENCES:
%[1] W. J. Ellison, "Permittivity of pure water, at standard atmospheric
%    pressure, over the frequency range 0-25 THz and the temperature range
%    0-100 °C," Journal of Physical and Chemical Reference Data, vol. 36,
%    no. 1, Mar. 2007.
%[2] T. Meissner and F. J. Wentz, "The complex dielectric constant of pure
%    water from microwave satellite observations," IEEE Transactions on
%    Geoscience and Remote Sensing, vol. 42, no. 9, pp. 1836-1849, Sep.
%    2004.
%[3] J. E. K. Laurens and K. E. Oughstun, "Electromagnetic impulse
%    response of triply-distilled water," in Proceedings of the 4th Ultra-
%    Wideband Short Pulse Electromagnetic Conference, Tel-Aviv, Israel,
%    14-19 Jun. 1999, pp. 243-264.
%[4] F. Shubitidze and U. Osterberg, "Phenomenological model to fit
%    complex permittivity data of water from radio to optical frequencies,"
%    Physical Review E, vol. 75, no. 4, Apr. 2007.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<1||isempty(T))
   T=15;%Centigrade
end

if(nargin<2||isempty(v))
    v=2e9;%2 GHz
end

if(nargin<3||isempty(algorithm))
    if(v<=25e12)
        algorithm=0;
    else
        algorithm=2;
    end
end

switch(algorithm)
    case 0
        e=EllisonSol(T,v);
    case 1
        e=MeissnerSol(T,v);
    case 2
        e=LaurensSol(v);
    case 3
        e=ShubitidzeSol(T,v);
    otherwise
        error('Unknown algorithm specified.')
end

end

function e=MeissnerSol(T,v)
%MEISSNERSOL This function implements the solution for the complex
%           permittivity of freshwater in [1], which is valid for
%           temperatures between -20 (supercooled water) and 40 degrees
%           centigrade for frequencies up to 500GHz.
%
%The sign of the complex part is flipped here, since in [1], they use a
%negative sign for the complex paart and in this function absorption losses
%are given a positive sign.
%
%REFERENCES
%[1] T. Meissner and F. J. Wentz,"The complex dielectric constant of pure
%    water from microwave satellite observations," IEEE Transactions on
%    Geoscience and Remote Sensing, vol. 42, no. 9, pp. 1836-1849, Sep.
%    2004.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.

%Convert from Hz to GHz.
v=v/1e9;

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

%Equation 8 in [1].
e1=a0+a1*T+a2*T^2;
v1=(45+T)/(a3+a4*T+a5*T^2);
eInf=a6+a7*T;
v2=(45+T)/(a8+a9*T+a10*T^2);

%Equation 7 in [1].
eS=(3.70886e4-8.2168e1*T)/(4.21854e2+T);

%Equation 6 in 1. The sigma term is 0 for freshwater.
e=(eS-e1)/(1+1i*v/v1)+(e1-eInf)/(1+1i*v/v2)+eInf;

end

function e=EllisonSol(t,v)
%%ELLISONSOL This function implements the solution for the complex
%            permittivity of freshwater in [1], which is valid for
%            temperatures between -5 and 100 degrees centigrade and
%            frequencies from 0 to 25THz.
%
%t is in degrees Centigrade and v is in Hz.
%
%REFERENCES:
%[1] W. J. Ellison, "Permittivity of pure water, at standard atmospheric
%    pressure, over the frequency range 0-25 thz and the temperature range
%    0-100 °C," Journal of Physical and Chemical Reference Data, vol. 36,
%    no. 1, Mar. 2007.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.

%The values below are taken from Table 2.
%Parameters for the three relaxation terms.
a1=79.23882;
a2=3.815866;
a3=1.634967;
tc=133.1383;
b1=0.004300598;
b2=0.01117295;
b3=0.006841548;
c1=1.382264e-13;
c2=3.510354e-16;
c3=6.30035e-15;
d1=652.7648;
d2=1249.533;
d3=405.5169;

%First resonance parameters
p0=0.8379692;
p1=-0.006118594;
p2=-0.000012936798;
p3=4235901000000;
p4=-14260880000;
p5=273815700;
p6=-1246943;
p7=9.618642e-14;
p8=1.795786e-16;
p9=-9.310017e-18;
p10=1.655473e-19;

%Second resonance parameters
p11=0.6165532;
p12=0.007238532;
p13=-0.00009523366;
p14=15983170000000;
p15=-74413570000;
p16=497448000;
p17=2.882476e-14;
p18=-3.142118e-16;
p19=3.528051e-18;

es=87.9144-0.404399*t+9.58726e-4*t^2-1.32802e-6*t^3;
delta1=a1*exp(-b1*t);
delta2=a2*exp(-b2*t);
delta3=a3*exp(-b3*t);
tau1=c1*exp(d1/(t+tc));
tau2=c2*exp(d2/(t+tc));
tau3=c3*exp(d3/(t+tc));

delta4=p0+p1*t+p2*t^2;
f0=p3+p4*t+p5*t^2+p6*t^3;
tau4=p7+p8*t+p9*t^2+p10*t^3;
delta5=p11+p12*t+p13*t^2;
f1=p14+p15*t+p16*t^2;
tau5=p17+p18*t+p19*t^2;

%Equation 17a
ep=es-(2*pi*v)^2*(tau1^2*delta1/(1+(2*pi*v*tau1)^2)...
    +tau2^2*delta2/(1+(2*pi*v*tau2)^2)...
    +tau3^2*delta3/(1+(2*pi*v*tau3)^2))...
    -(2*pi*tau4)^2*(delta4/2)*(v*(f0+v)/(1+(2*pi*tau4*(f0+v))^2)...
    -v*(f0-v)/(1+(2*pi*tau4*(f0-v))^2))...
    -(2*pi*tau5)^2*(delta5/2)*(v*(f1+v)/(1+(2*pi*tau5*(f1+v))^2)...
    -v*(f1-v)/(1+(2*pi*tau5*(f1-v))^2));

%Equaion 17b
epp=2*pi*v*(tau1*delta1/(1+(2*pi*v*tau1)^2)...
    +tau2*delta2/(1+(2*pi*v*tau2)^2)...
    +tau3*delta3/(1+(2*pi*v*tau3)^2))...
    +pi*v*tau4*delta4*(1/(1+(2*pi*tau4*(f0+v))^2)...
    +1/(1+(2*pi*tau4*(f0-v))^2))...
    +pi*v*tau5*delta5*(1/(1+(2*pi*tau5*(f1+v))^2)...
    +1/(1+(2*pi*tau5*(f1-v)^2)));


%The sign of the imaginary part is flipped compared to the paper.
e=ep-1i*epp;

end

function e=LaurensSol(f)
%%LAURENSSOL This function implements the solution for the complex
%            permittivity of freshwater in [1]. The valid temperature range
%            is not specified in [1]. The valid frequency range seems to
%            span 1e10 to 1e18Hz.
%
%REFERENCES:
%[1] J. E. K. Laurens and K. E. Oughstun, "Electromagnetic impulse
%    response of triply-distilled water," in Proceedings of the 4th Ultra-
%    Wideband Short Pulse Electromagnetic Conference, Tel-Aviv, Israel,
%    14-19 Jun. 1999, pp. 243-264.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.

twoPi=2*pi;

omega=f*twoPi;

%Values from table 1
eInf=1;
a0=74.1;
a2=2.90;
tau0=8.44e-12;
tau2=6.05e-14;
tauf0=4.93e-14;
tauf2=8.59e-15;
omega11=1.8e13*twoPi;
omega13=4.9e13*twoPi;
omega15=1e14*twoPi;
omega17=3.7e15*twoPi;
b11=1.2e13*twoPi;
b13=6.8e12*twoPi;
b15=2.0e13*twoPi;
b17=3.2e15*twoPi;
delta11=4.3e12*twoPi;
delta13=8.4e11*twoPi;
delta15=2.8e12*twoPi;
delta17=8.0e14*twoPi;

%Equation 8 in [1].
e=eInf+a0/((1-1i*omega*tau0)*(1-1i*omega*tauf0))...
    +a2/((1-1i*omega*tau2)*(1-1i*omega*tauf2))...
    -b11^2/(omega^2-omega11^2+2*1i*delta11*omega)...
-b13^2/(omega^2-omega13^2+2*1i*delta13*omega)...
-b15^2/(omega^2-omega15^2+2*1i*delta15*omega)...
-b17^2/(omega^2-omega17^2+2*1i*delta17*omega);

%Flip the sign on the imaginary part compared to the paper.
e=real(e)-1i*imag(e);

end

function e=ShubitidzeSol(T,f)
%%SHUBITIDZESOL This implements the solution in [1] using their data and
%       the specific constants that they used. The data is not available in
%       the paper and was obtained by email from the authors. Some of the
%       equations that they used differ from the paper and are commented.
%
%REFERENCES:
%[1] F. Shubitidze and U. Osterberg, "Phenomenological model to fit
%    complex permittivity data of water from radio to optical frequencies,"
%    Physical Review E, vol. 75, no. 4, Apr. 2007.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.

%For some reason they fit everything to the frequency divided by 3 rather
%than the actual frequency.
f=f/3;

%Convert the temperature from Centigrade to Kelvin.
T=273.15+T;%The temperature in degrees Kelvin.

Dn=[6; 6; 6; 8; 8;12;12;10; 8;-8; 8;12;12;12;12];
Nn=[5; 5; 5; 3; 3; 3; 3; 3; 3;-3; 3; 3; 5; 2; 4];

%Allocate space for the tables
D=NaN(15,12);
N=NaN(15,5);
    
D(1, 1:6)= [ 0.476809954333578E-08-1i*0.173232046937905E-17,...
            -0.719099595872712E-06+1i*0.122081819893495E-15,...
             0.477692820332839E-04-1i*0.137524195542974E-14,...
            -0.163009397116141E-02-1i*0.475945490194174E-13,...
             0.296319139300348E-01+1i*0.107802192241306E-11,...
            -0.272411400873583E+00-1i*0.621655670276182E-11];
D(2, 1:6)= [ 0.241131336175290E-04-1i*0.164268975137902E-14,...
            -0.836903227940087E-03+1i*0.970762542481097E-14,...
             0.122375268871338E-01+1i*0.259036024827046E-12,...
            -0.961697606072290E-01-1i*0.327596108637332E-11,...
             0.426987493310241E+00+1i*0.135950654483177E-10,...
            -0.101266213409972E+01-1i*0.193957599129896E-10];
D(3, 1:6)= [ 0.130017472062263E-02+1i*0.896125240401395E-14,...
            -0.236305125856654E-01-1i*0.134586302375215E-12,...
             0.178892671097471E+00+1i*0.808218895777117E-12,...
            -0.722071179402915E+00-1i*0.242590958103592E-11,...
             0.163896035202330E+01+1i*0.363953905923064E-11,...
            -0.198354687663382E+01-1i*0.218343743306747E-11];
D(4, 1:8)= [ 0.155735966726143E-01+1i*0.212617232756609E-11,...
            -0.218232732023098E+00-1i*0.238102244371992E-10,...
             0.132274747175931E+01+1i*0.112370398234050E-09,...
            -0.452882728330265E+01-1i*0.289769844491498E-09,...
             0.957932878800585E+01+1i*0.440944443584403E-09,...
            -0.128181082133538E+02-1i*0.395885288313788E-09,...
             0.105969227296838E+02+1i*0.194122288709777E-09,...
            -0.494924884090296E+01-1i*0.400965673243015E-10];
D(5, 1:8)= [ 0.202460765781801E+02+1i*0.535186835448412E-08,...
            -0.114183858317919E+03-1i*0.253929363072498E-07,...
             0.279606928228480E+03+1i*0.512629329092922E-07,...
            -0.388286960790323E+03-1i*0.570939401885430E-07,...
             0.334462529770788E+03+1i*0.378919439430284E-07,...
            -0.183000258194299E+03-1i*0.149856626963725E-07,...
             0.621145890179884E+02+1i*0.326985131023317E-08,...
            -0.119588738487024E+02-1i*0.303664602717245E-09];
D(6, 1:12)=[-0.108878585945448E+07+1i*0.139529046502617E-03,...
             0.165401894507161E+04-1i*0.254728044905815E-04,...
             0.491000361056760E+07-1i*0.465818254481910E-03,...
            -0.852263407934120E+07+1i*0.749032697650696E-03,...
             0.746324485698593E+07-1i*0.581882248965049E-03,...
            -0.411692810485142E+07+1i*0.276816322398357E-03,...
             0.154704195985514E+07-1i*0.871105313272165E-04,...
            -0.408576204542642E+06+1i*0.185890724589496E-04,...
             0.762227429601001E+05-1i*0.267041614636237E-05,...
            -0.986994440748875E+04+1i*0.247758154739893E-06,...
             0.846128929909977E+03-1i*0.134187114889155E-07,...
            -0.432614383723785E+02+1i*0.322173542013347E-09];
D(7, 1:12)=[ 0.803005744541716E+09+1i*0.121377235016101E+00,...
            -0.656169794501136E+04-1i*0.170705164288651E-02,...
             0.173526744337001E+07+1i*0.934001079097131E-01,...
            -0.403324409742438E+09-1i*0.150835526011787E+00,...
             0.354957921050070E+09+1i*0.886570277056129E-01,...
            -0.147833304290869E+09-1i*0.285232247781373E-01,...
             0.373113801747888E+08+1i*0.572571491380883E-02,...
            -0.620892733735174E+07-1i*0.753102083919981E-03,...
             0.701129374961377E+06+1i*0.652800115013268E-04,...
            -0.534654524138770E+05-1i*0.360694458661981E-05,...
             0.264614121429424E+04+1i*0.115479854079454E-06,...
            -0.769490483116168E+02-1i*0.163403257806549E-08];
D(8, 1:10)=[-0.894578315069949E+16-1i*0.204116377664408E+07,...
            -0.228039873886267E+11-1i*0.345985653273340E+01,...
             0.211988152149394E+15+1i*0.139826580977095E+05,...
            -0.250464654847871E+14-1i*0.737352998830355E+03,...
             0.140236215802872E+13+1i*0.212892609702946E+02,...
            -0.463557565191333E+11-1i*0.390585363145073E+00,...
             0.966305113030240E+09+1i*0.456272265440882E-02,...
            -0.128646157525605E+08-1i*0.329648621158009E-04,...
             0.106214487648793E+06+1i*0.134217358382534E-06,...
            -0.495751414292685E+03-1i*0.235610588323351E-09];
D(9, 1:8)= [ 0.149019866683351E-17-1i*0.122271135472144E-27,...
            -0.124957795120504E-14-1i*0.982136497307247E-25,...
             0.229172602244895E-11-1i*0.800241262779948E-22,...
            -0.838253461819030E-09+1i*0.103351991558314E-19,...
             0.183290973858022E-06-1i*0.121352834043285E-18,...
            -0.237343732278274E-04+1i*0.160710666703582E-17,...
             0.173476115067901E-02-1i*0.300174170921429E-14,...
            -0.656704860904070E-01+1i*0.146060356989103E-12];
D(10, 1:8)=[ 0.178809247227920E-17-1i*0.175284368007891E-27,...
            -0.267255639603626E-14-1i*0.170286304219136E-24,...
             0.312869505635337E-11-1i*0.105785044354370E-21,...
            -0.118495273918795E-08+1i*0.648534483140379E-20,...
             0.247960248087857E-06+1i*0.118249690754574E-17,...
            -0.297217070827536E-04-1i*0.114446435744309E-15,...
             0.200503231992348E-02+1i*0.124279100112695E-14,...
            -0.704384451157339E-01+1i*0.804176399041238E-13];
D(11, 1:8)=[ 0.204912166854986E-09+1i*0.740946483161539E-21,...
            -0.300723722967368E-07+1i*0.121225361345385E-18,...
             0.185939189671353E-05-1i*0.167487214880736E-16,...
            -0.632461441341890E-04+1i*0.715518183086192E-15,...
             0.129760826259057E-02-1i*0.149710415180957E-13,...
            -0.164872427141606E-01+1i*0.167085292100706E-12,...
             0.127071217403121E+00-1i*0.949953869968032E-12,...
            -0.544988200106952E+00+1i*0.213641761478217E-11];
D(12, 1:12)=[0.299628414550434E+10-1i*0.334110360137443E-01,...
             0.145690563402376E+07-1i*0.148975086720719E-04,...
             0.710909694842852E+09+1i*0.220421605989130E-01,...
            -0.184063965144037E+10-1i*0.144383852876144E-01,...
             0.121371874606673E+10+1i*0.472215851784406E-02,...
            -0.417870056901506E+09-1i*0.959668372370581E-03,...
             0.895729407618336E+08+1i*0.130440605503809E-03,...
            -0.127965692123849E+08-1i*0.121443867271552E-04,...
             0.124679945840470E+07+1i*0.768767506880346E-06,...
            -0.822377579508378E+05-1i*0.317564857190269E-07,...
             0.352476589406243E+04+1i*0.774324213282137E-09,...
            -0.888073935742613E+02-1i*0.847262890687834E-11];
D(13, 1:12)=[0.224199656690788E-06-1i*0.122076571674529E-16,...
            -0.102964281021972E-04+1i*0.486228296276323E-15,...
             0.214080281848989E-03-1i*0.867919275863108E-14,...
            -0.266445928043494E-02+1i*0.916004042487220E-13,...
             0.221095210138599E-01-1i*0.634742608372072E-12,...
            -0.128872922172256E+00+1i*0.303073109904809E-11,...
             0.541142844180795E+00-1i*0.101711590989598E-10,...
            -0.164966299011010E+01+1i*0.239878141188707E-10,...
             0.362435984839570E+01-1i*0.389597603979899E-10,...
            -0.559809399336247E+01+1i*0.415031002544452E-10,...
             0.577171759884929E+01-1i*0.261020665612764E-10,...
            -0.356747456162304E+01+1i*0.734320582787060E-11];
D(14, 1:12)=[0.341093568369667E+04-1i*0.510460492874808E-03,...
            -0.239828324873696E+05+1i*0.114166183344234E-02,...
             0.760675085649523E+05-1i*0.202039291374953E-02,...
            -0.140672773341506E+06+1i*0.260118786129950E-02,...
             0.170839532080944E+06-1i*0.225261076196367E-02,...
            -0.143551845627351E+06+1i*0.138730917218846E-02,...
             0.856203252973933E+05-1i*0.625423008479579E-03,...
            -0.365563664270621E+05+1i*0.206975294778838E-03,...
             0.111014926075174E+05-1i*0.490972809165942E-04,...
            -0.234146538646201E+04+1i*0.788549298186779E-05,...
             0.325993884152267E+03-1i*0.764282384936678E-06,...
            -0.269347765388244E+02+1i*0.335667647269352E-07];
D(15, 1:12)=[0.238335596359778E+07-1i*0.133154831632230E+00,...
             0.391449227591339E+04+1i*0.214205025642204E-03,...
             0.837374868764651E+06+1i*0.324001933626814E-01,...
            -0.661443215505139E+07+1i*0.144912904620094E+00,...
             0.851291556047272E+07-1i*0.207927740986336E+00,...
            -0.550435667011906E+07+1i*0.126616847570447E+00,...
             0.219368084346928E+07-1i*0.448098242819989E-01,...
            -0.580873356841553E+06+1i*0.100828044082053E-01,...
             0.104841358043345E+06-1i*0.147214587241373E-02,...
            -0.128172737596404E+05+1i*0.135781993842008E-03,...
             0.101934594834842E+04-1i*0.721584496852315E-05,...
            -0.477210963496044E+02+1i*0.168848958419613E-06];
N(1, 1:5)= [ 0.182974540514735E-17+1i*0.436532247091081E-08,...
             0.492148609773319E-15+1i*0.807839120418427E-07,...
            -0.460976537790477E-13-1i*0.960987981913376E-05,...
             0.108052354980312E-11+1i*0.152531943801301E-03,...
            -0.673990376677299E-11+1i*0.815016650912784E-05];
N(1, 1:5)=0.5*N(1, 1:5);
N(2, 1:5)= [ 0.141642275416902E-13+1i*0.117849098290432E-04,...
            -0.298662873436138E-12-1i*0.256993512248829E-03,...
             0.235195240084723E-11+1i*0.209849712989052E-02,...
            -0.819755381538116E-11-1i*0.760599149713317E-02,...
             0.106714990884535E-10+1i*0.103299148881628E-01];
N(3, 1:5)= [ 0.225776302560162E-13+1i*0.297597838715721E-05,...
            -0.280698824713530E-12-1i*0.493018161633052E-04,...
             0.130530641713520E-11+1i*0.313377301298601E-03,...
            -0.269087427816112E-11-1i*0.876523213775654E-03,...
             0.207502384573184E-11+1i*0.897778616453127E-03];
N(4, 1:3)= [ 0.276732348106894E-14+1i*0.365818723726163E-06,...
            -0.951405243058516E-14-1i*0.137179627633014E-05,...
             0.821058224250961E-14+1i*0.138190735082626E-05];
N(5, 1:3)= [ 0.127343414569749E-11+1i*0.110150555057169E-03,...
            -0.190745404345165E-11-1i*0.165691205763030E-03,...
             0.717626083163194E-12+1i*0.631681186214675E-04];
N(6, 1:3)= [-0.466034743567454E-05+1i*0.470185088300228E+03,...
             0.285210102956482E-05-1i*0.396482545039419E+03,...
            -0.372812008661821E-06+1i*0.578527956264391E+02];
N(7, 1:3)= [-0.187536896786991E-04+1i*0.215376481837741E+03,...
             0.555747591734941E-05-1i*0.101365402220132E+03,...
            -0.407824254637021E-06+1i*0.112457765722161E+02];
N(8, 1:3)= [ 0.265717800281251E+07-1i*0.104445524697657E+16,...
            -0.171144144134818E+06+1i*0.753419010249941E+14,...
             0.171419528179475E+04-1i*0.265261899658379E+13];
N(9, 1:3)= [-0.235790443468372E-19+1i*0.109002749580609E-28,...
             0.696667500461018E-15-1i*0.115779728104005E-24,...
            -0.313729673946882E-13+1i*0.921628159981826E-23];
N(10, 1:3)=[ 0.231485323666048E-19+1i*0.162830578428740E-28,...
             0.463538768577934E-15-1i*0.173860802216102E-24,...
            -0.226738613227731E-13+1i*0.998169043270983E-23];
N(11, 1:3)=[-0.276390114210373E-12+1i*0.618241519954187E-22,...
             0.148091770152658E-10-1i*0.341791265290408E-20,...
            -0.182538556631705E-09+1i*0.471709115833368E-19];
N(12, 1:3)=[ 0.272214752606561E+05+1i*0.148991756260095E-05,...
            -0.862450002489323E+04-1i*0.411780972780551E-06,...
             0.669203005970406E+03+1i*0.284544146142263E-07];
N(13, 1:5)=[ 0.153533665065256E-10-1i*0.104634484118712E-19,...
            -0.264256220921621E-09+1i*0.190269669907308E-18,...
             0.164389800782089E-08-1i*0.123901837111523E-17,...
            -0.435961876036075E-08+1i*0.340326127817120E-17,...
             0.417910289059908E-08-1i*0.334814804589455E-17];
N(14, 1:2)=[-0.173516177369559E+02-1i*0.980122702424710E-04,...
             0.927970308968244E+02+1i*0.614408926035496E-04];
N(15, 1:4)=[ 0.198345596413140E+04+1i*0.510830645059674E-03,...
            -0.105991917521582E+04-1i*0.268357297744296E-03,...
             0.188830571272145E+03+1i*0.469692453990141E-04,...
            -0.112155666428697E+02-1i*0.273885087933284E-05];
        
%They use the frequency, not the frequency multiplied by 2*pi, for most
%things.
omega=f;
%The above coefficients below were fit to a scaled frequency, not
%the actual frequency in Hertz.
f=f/1e14;

%Evaluate Equation 17.
aOmega=0;
for i=1:15
    %The sum in the denominator of Equation 17.
    DSum=f.^Dn(i);
    for id=0:abs(Dn(i))-1
        DSum=DSum+D(i,id+1)*(f.^id);
    end

    %The sum in the numerator of Equation 17.
    NSum=0;
    for in=0:abs(Nn(i))-1
        NSum=NSum+N(i,in+1)*(f.^in);
    end
    %The a term in Equation 17.
    aOmega=aOmega+NSum./DSum;
end

h=6.62*1e-34;%Approximate value of Planck's constant.
k=1.38e-23;%Approximate value of Boltzman's constant.
%The A and B constants in Equations 15 and 19.
A_u=0.37122e-9;
B_u=0.3628e-9;
Me=9.31e-31;%Approximate value of the mass of an electron.

%From Equation 16 in the paper. This is an instance where omega is used as
%an angular frequency. For the rest of their fiting, they took omega to be
%the actual frequency in Hertz.
r=sqrt(h/Me)./sqrt(omega*2*pi^2);

%Ep from Equation 15.
Ep=((B_u./(r)).^12-(A_u./(r)).^6);
%Two times Equation 13, because that is what they fit to, not to Equation
%13.
F=2./(1+exp(h*omega/(1.5*k*T)-Ep));

%The omega_0 constant in Equation 12.
omega0=3.0e+15;
%The gamma damping constants in Equation 12.
gamma=2*3*1.198e+15;
%The es constant for the Lorentz distribution in Equation 12
es=2.5;            
%The Lorentz distribution's high frequency response constant in Equation
%12.
eInf=1;

%The high-frequency response constant in Equation 11 (in the Debye
%distribution).
Xi0=7.4;
%The two tau constants in Equation 11.
tau1=0.25*1e-10;
tau2=3*0.1e-9;

%Equation 11. The weird floor operation in Equation 11 in the paper
%shouldn't be there. Also, the logarithm of the tau ratio should only be a
%single logarithm, not a double one.
Debye=Xi0*(1-1/log(tau2/tau1)*log((1i*omega*tau2+1)./(1i*omega*tau1+1))).*sqrt(F);
%Equation 12, except the omega^2 term is multiplied by 9, because that is
%how they actually fit it.
Lorentz=sqrt(eInf-omega0^2*(es-eInf)./(omega.^2*9-2*1i*gamma*omega-omega0^2).*F);
%Equation 10
bOmega=Debye+Lorentz;

%Equation 3, as they actually implemented it, not as it is in the paper.
eReal=(1+real(aOmega)).*real(bOmega);
%Equation 4.
eImag=imag(aOmega).*imag(bOmega);

%The paper computes a complex index of refraction. We want the complex
%relative permittivity. Approximating the permeability as 1, we just need
%to square the index of refraction to get the permittivity.
e=(eReal+1i*eImag).^2;

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
