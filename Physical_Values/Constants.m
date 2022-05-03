classdef Constants
%%CONSTANTS Static methods for obtaining standard physical constants.
%
%Constants can be accessed using Constants.constantName without
%instantiating this class.
%
%The WGS84 properties are from  [1]. The formula for the mean radius is
%from [14] and is also used with the WGS72 and GRS80 proeprties.
%
%The WGS72 properties are from [2].
%
%The GRS80 properties are from [6].
%
%Tha IERS' mean Earth rotation rate, which plays a role in computations of
%the length-of-day Earth orientation parameter is from
%http://hpiers.obspm.fr/eop-pc/earthor/ut1lod/UT1.html
%
%The value of the amount subjected from a Julian date to get a Modified
%Julian date is given in the IERS conventions [6]. The number of seconds in
%a TT Julian day is given as 86400 in [6], indicated as the astronomical
%unit of time.
%
%The ellipsoid properties for the EGM2008 model are taken from 
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/README_FIRST.pdf
%where the value of GM and the semi-major axis of the reference ellipsoid
%differ from the WGS84 values. Note that these values are the same for the
%EGM96 gravitational model.
%
%The value of the International Geomagnetic Reference Field's (IGRF's)
%reference sphere radius, which is used in the 11th edition of the IGRF, is
%from [3].
%
%The value of the reference sphere radius in the World Magnetic Model for
%the year 2010 is given in [4].
%
%The CODATA 2014 recommended values are from [5]. The data is also on the
%site http://physics.nist.gov/cuu/Constants/ .
%
%Various physical parameters for space and the Earth are from [6]. This
%includes the drift rates for the C20, C30 and C40 terms in the EGM2008
%model, which are given in table 6.2 on page 80.
%
%The values for the gravitational parameters of the Moon are taken from the
%GL0900C, which can be obtained from
%http://pds-geosciences.wustl.edu/missions/grail/default.htm
%and which is very briefly documented in [7] and which is consistent with
%the DE430 ephemerides.
%
%The parameters for the JPL reference ellipsoid of the moon are taken from
%[8].
%
%The SI standard values are from the 8th Edition of the International
%System of Units at
%http://www.bipm.org/en/si/si_brochure/
%
%ICAO's International Standard Atmosphere is used for the standard
%temperature, pressure and relative humidity.
%
%The atomic weights are from the 2013 Table of Standard Atomic Weights that
%is published by the Commission on Isotopic Abundances and Atomic Weights,
%Commission II.I of the International Union of Pure and Applied Chemistry
%and is available at
%http://www.ciaaw.org/publications.htm
%The atomic weights are given for constitutents of gasses that tend to be
%in the atmosphere. Sometimes, the CIAAW gives a range of values for an
%element instead of just one value. In such an instance, the endpoints of
%the range are just averaged here.
%
%The method gasProp returns the atomic weight for a given substance. If
%temperatures are provided, then the ideal gas specific heat and the second
%virial coefficient and its first two derivatives for the specified
%substance are returned as well. The data for all of the second virial
%coefficients except Helium and water were was taken from [9]
%
%The second virial coefficient for Helium is taken from cubic spline
%interpolation of the tabulated values in [10]. The second virial
%coefficient for water was taken from [11] where it was assumed that the
%valid range matched the temperature range studied in the paper.
%The data for all of the ideal gas specific heats at a constant pressure
%were taken from [12].
%
%The Euler-Mascheroni Constant is as mentioned in [13]. Khinchin's constant
%is described in [14]. The golden ratio and golden angles are described in
%[16] and [17].
%
%REFERENCES:
%[1] Department of Defense, "Department of Defense world geodetic system
%    1984: Its definition and relationships with local geodetic systems,"
%    National Imagery and Mapping Agency, Tech. Rep. NIMA TR8350.2, Jun.
%    2004, third Edition, Amendment 2. [Online]. Available:
%    http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
%[2] World Geodetic System Committee, "The Department of Defense world
%    geodetic system 1972," Defense Mapping Agency, Washington, D.C., Tech.
%    Rep., May 1974. [Online]. Available:
%    http://www.dtic.mil/dtic/tr/fulltext/u2/a110165.pdf
%[3] International Association of Geomagnetism and Aeronomy, Working
%    Group V-MOD, "International geomagnetic reference field: Eleventh
%    generation," Geophysical Journal International, vol. 183, no. 3, pp.
%    1216-1230, Dec. 2010.
%[4] S. Maus, S. McLean, M. Nair, and C. Rollins, "The US/UK
%    world magnetic model for 2010-2015," National Oceanographic and
%    Atmospheric Organization, Tech. Rep. NESDIS/NGDC, 2010. [Online].
%    Available: http://www.ngdc.noaa.gov/geomag/WMM/
%[5] P. J. Mohr, D. B. Newell, and B. N. Taylor, "CODATA recommended
%    values of the fundamental physical constants: 2014," ArXiv
%    21 Jul. 2015. [Online]. Available: http://arxiv.org/abs/1507.07956
%[6] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%[7] A. S. Konopliv, R. S. Park, D.-N. Yuan, S. W. Asmar, and et. al.,
%    "High-resolution lunar gravity fields from the GRAIL primary and
%    extended missions," Geophysical Research Letters, vol. 41, no. 5, pp.
%    1452-1458, 16 Mar. 2014.
%[8] R. B. Roncoli, "Lunar constants and models document," Jet Propulsion
%    Laboratory, California Institute of Technology, Tech. Rep. JPL D-32296,
%    23 Sep. 2005. [Online].
%    Available: http: //www.hq.nasa.gov/alsj/lunar cmd 2005 jpl d32296.pdf
%[9] D. Ambrose, M. B. Ewing, and M. L. McGlashan. (2014, Mar.) Kaye & Laby
%    tables of physical & chemical constants. National Physical Laboratory.
%    [Online]. Available: http://www.kayelaby.npl.co.uk/chemistry/3_5/3_5.html
%[10] J. M. H. Levelt Sengers, M. Klein, and J. S.Gallagher, "Pressure-
%    Volume-Temperature Relationships of Gases; Virial Coefficients," in
%    American Institute of Physics Handbook, edited by D. E. Grey (McGraw-
%    Hill, New York, 1972), 3rd ed., Chap. 4i, pp. 4-204-4-221.
%[11] R. W. Hyland, "A correlation for the second interaction virial
%    coefficients and enhancement factors for moist air," Journal of
%    Research of National Bureau of Standards - A Physics and Chemistry,
%    vol. 79A, no. 4, pp. 551-560, Jul. - Aug. 1975.
%[12] Y. S. Touloukian and T. Makita, "Specific heat nonmetallic liquids
%    and gasses," in Thermophysical Properties of Matter. New York:
%    IFI/Plenum, 1970, vol. 6.
%[13] Weisstein, Eric W. "Euler-Mascheroni Constant." From MathWorld--A
%    Wolfram Web Resource. http://mathworld.wolfram.com/Euler-MascheroniConstant.html
%[14] H. Mortiz, "Geodetic Reference System 1980," Bulletin G�od�sique,
%    vol. 54, no. 3, pp. 395-405, Sep. 1980. Given with corrections at
%    https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf
%[15] Weisstein, Eric W. "Khinchin's Constant." From MathWorld--A 
%     Wolfram Web Resource. https://mathworld.wolfram.com/KhinchinsConstant.html
%[16] Weisstein, Eric W. "Golden Ratio." From MathWorld--A Wolfram Web
%     Resource. https://mathworld.wolfram.com/GoldenRatio.html
%[17] Weisstein, Eric W. "Golden Angle." From MathWorld--A Wolfram Web
%     Resource. https://mathworld.wolfram.com/GoldenAngle.html
%
%January 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

properties (Constant)
    %WGS84 Properties
    %GM is the universal gravitation constant times the mass of the Earth.
    WGS84GMWithAtmosphere=3.986004418*10^(14);%m^3/s^2
    WGS84GMWithoutAtmosphere=3.9860009*10^(14);%m^3/s^2
    WGS84EarthRotationRate=7292115.0*10^(-11);%radians per second.
    %The following 4 parameters should be consistent with values in the
    %IAU's SOFA library.
    WGS84SemiMajorAxis=6378137.00;%m
    WGS84InverseFlattening=298.257223563;%Unitless
    WGS84Flattening=1/Constants.WGS84InverseFlattening;%Unitless
    WGS84SemiMinorAxis=(Constants.WGS84SemiMajorAxis*(Constants.WGS84InverseFlattening-1))/Constants.WGS84InverseFlattening;
    %The squared first numerical eccentricity of the ellipsoid.
    WGS84e2=2*Constants.WGS84Flattening-Constants.WGS84Flattening^2;
    
    %This is consistent with the value in Table 3.3 of [1]. The formula to
    %derive the mean radius is given in [14].
    WGS84MeanRadius=(2*Constants.WGS84SemiMajorAxis+Constants.WGS84SemiMinorAxis)/3;
    
    %WGS72 Properties
    WGS72GMWithAtmosphere=398600.8*10^(9);%m^3/s^2
    WGS72GMWithoutAtmosphere=398600.5*10^(9);%m^3/s^2
    WGS72EarthRotationRate=0.7292115147*10^(-4);%radians per second.
    WGS72SemiMajorAxis=6378135;%m
    WGS72InverseFlattening=298.26;%Unitless
    WGS72Flattening=1/Constants.WGS72InverseFlattening;%Unitless
    WGS72C20Bar=-484.1605e-6;%Unitless, tide free
    WGS72SemiMinorAxis=(Constants.WGS72SemiMajorAxis*(Constants.WGS72InverseFlattening-1))/Constants.WGS72InverseFlattening;
    WGS72MeanRadius=(2*Constants.WGS72SemiMajorAxis+Constants.WGS72SemiMinorAxis)/3;
    %The squared first numerical eccentricity of the ellipsoid.
    WGS72e2=2*Constants.WGS72Flattening-Constants.WGS72Flattening^2;
    
    %GRS80 Properties
    GRS80GMWithAtmosphere=3.986005e14;%m^3/s^s
    GRS80SemiMajorAxis=6378137;%m
    GRS80InverseFlattening=298.257222101;%Unitless
    GRS80Flattening=1/Constants.GRS80InverseFlattening
    GRS80J2=1.08263e-3;%Dynamical form factor
    GRS80EarthRotationRate=7.292115e-5;%radians/second
    GRS80SemiMinorAxis=(Constants.GRS80SemiMajorAxis*(Constants.GRS80InverseFlattening-1))/Constants.GRS80InverseFlattening;
    GRS80MeanRadius=(2*Constants.GRS80SemiMajorAxis+Constants.GRS80SemiMinorAxis)/3;
    %The squared first numerical eccentricity of the ellipsoid.
    GRS80e2=2*Constants.GRS80Flattening-Constants.GRS80Flattening^2;
    
    %This is the mean rotation rate at epoch 1820, from which the
    %Length-of-day Earth orientation parameter is related.
    IERSMeanEarthRotationRate=72921151.467064e-12;%Radians per second
    
    %This is the amount subtracted from a Julian date to transform it into
    %a Modified Julian date. from 
    MJDOffset=2400000.5;
    %TT,TAI Julian days have no leap seconds.
    secondsPerTTJulianDay=86400;

    %The EGM2008 model; the same values are used for the EGM96 model. These
    %values are needed for using the spherical harmonic coefficients in
    %the models. They are also the defining parameters of the reference
    %ellipsoid used in the model
    EGM2008GM=3986004.415*10^8;%m^3/s^2
    EGM2008SemiMajorAxis=6378136.3%m
    EGM2008EarthRotationRate=7292115*10^(-11);%rad/s
    EGM2008C20BarRefEllips=-484.1654767*10^(-6);%Unitless, tide-free,
    %defines the reference ellipsoid.
    
    %The following are drift rates in the fully normalized spherical
    %harmonic constants in the EGM2008 model. The units is the change per
    %year (since the spherical harmonic coefficients themselves are
    %unitless). The epoch date for the EGM2008 model is J2000.0. The change
    %is assumed to be in Julian years terrestrial time (365.25 days per
    %year with a day having exactly 86400 seconds); however, the precision
    %of the model is low enough that UTC years and fractional years
    %can be used and that leap seconds should not matter.
    EGM2008C20BarDot=11.6e-12;
    EGM2008C30BarDot=1.3e-11/sqrt(2*3+1);
    EGM2008C40BarDot=1.4e-11/sqrt(2*4+1);
    
    EGM96GM=3986004.415*10^8;%m^3/s^2
    EGM96SemiMajorAxis=6378136.3%m
    EGM96EarthRotationRate=7292115*10^(-11);%rad/s
    
    %In the EGM96 model, coefficient rates per year are given for the fully
    %normalized forms of C20, C21 and S21 according to
    %http://cddisa.gsfc.nasa.gov/926/egm96/egm96.html
    %and in the readme file accompanying the EGM96 model.
    %These rates appear to have been abandoned in the EGM2008 model.
    EGM96C20BarDot=1.16275534e-11;%Per year since 1 January 1986
    EGM96C21BarDot=-0.32e-11;%Per year since 1 January 1986
    EGM96S21BarDot=1.62e-11;%Per year since 1 January 1986
    
    %For the EGM2008 model, when computing the gravitational potential of
    %the geoid, a different semi-major axis and flattening factor are used
    %than are used in the rest of the model. Also, the WGS-84 value of GM
    %is used along with the WGS-84 value of the Earth's rotation rate. This
    %is implied from
    %http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
    EGM2008GeoidFlattening=1/298.257686;%Unitless
    EGM2008GeoidSemiMajorAxis=6378136.58;%m
    
    %The universal gravitational constant in NASA JPL's GL0900C model,
    %converted from km^3/2^2 to m^3/s^2
    GL0900CGMMoon=0.49028001781922e+4*1000^3;
    %The radius of the Moon in NASA JPL's GL0900C model, converted from
    %kilometers to meters.
    GL0900CSemiMajorAxisMoon=0.1738e+4*1000;
    
    %These are parameters for a reference ellipsoid for the Moon.
    JPLMoonSemiMajorAxis=1737.4e3;%meters
    JPLMoonFlattening=0;
    %degrees per Julian TT day -> radians per second. --It is assumed that
    %"day" means Julian day TT.
    JPLMoonRotationRate=13.17635815*(pi/180)*(1/86400);
    
    %The radius of the IGRF11 reference sphere
    IGRF11SphereRad=6371.2*10^3;%meters
    
    %The radius of the reference sphere used in the WMM2010.
    WMM2010SphereRad=6371200;%meters
    
    %CODATA 2014 Recommended Values of the Fundamental Physical Constants
    speedOfLight=299792458;%(c) m/s
    magneticConstant=4*pi*10^(-7);%(mu_0) N/A^2
    electricConstant=1/(Constants.magneticConstant*Constants.speedOfLight^2);%(epsilon_0) F/m
    GravitationalConstant=6.67408e-11;%m^3/(kg*s^2)
    PlanckConstant=6.626070040e-34;%(h) Jule-seconds
    %The following is an alternative form of Planck's constant
    hBar=Constants.PlanckConstant/(2*pi);
    elementaryCharge=1.602176565e-19;%(e) Coulombs
    quantumMagneticFlux=Constants.PlanckConstant/(2*Constants.elementaryCharge);%Wb
    quantumConductance=2*Constants.elementaryCharge^2/Constants.PlanckConstant;%S
    electronMass=9.10938356e-31%(m_e) kilograms
    protonMass=1.672621898e-27;%(m_p) kilograms
    fineStructureConst=Constants.elementaryCharge^2/(4*pi*Constants.electricConstant*Constants.hBar*Constants.speedOfLight);%(alpha) Unitless.
    RydbergConstant=Constants.fineStructureConst^2*Constants.electronMass*Constants.speedOfLight/(2*Constants.PlanckConstant);%(R_Inf) meters
    AvogadroConstant=6.022140857e23;%(N_A) 1/mol
    FaradayConstant=Constants.AvogadroConstant*Constants.elementaryCharge;%(F) C/mol
    molarGasConstant=8.3144598;%(R) J/(mol K) (ideal gas constant)
    %The mean molecular mass of a substance in atomic mass units is the
    %mass in grams of one mole of the substance.
    BoltzmannConstant=Constants.molarGasConstant/Constants.AvogadroConstant;%(k) J/K
    StefanBoltzmannConstant=(pi^2/60)*Constants.BoltzmannConstant^4/(Constants.hBar^3*Constants.speedOfLight^2);%(sigma) W/(m^2*K^4)
    electronVolt=1.6021766208e-19;%(eV) Joules
    %Unified atomic mass unit
    atomicMassUnit=1.660539040e-27;%(u) kilograms
    
    %The mass ratios provided in ther IERS 2010 conventions.
    %From Table 1.1 of the IERS 2010 conventions.
    MoonEarthMassRatio=0.0123000371;
    %From Table 3.1 of the IERS 2010 conventions.
    SunEarthMassRatio=332946.048166;
    %Table 1.1 of the IERS 2010 conventions defines the equatorial radius
    %of the Earth.
    EarthEqRadius=6378136.6;%meters
    %From Table 1.1 of the IERS 2010 conventions.
    GMEarth=3.986004418*10^14;%m^3/s^2 Geocentric gravitational constant.
    %The universal gravitational constant according to the IERS 2010
    %conventions Table 1.1
    GravitationalConstantIERS=6.67428*10^(-11);%m^3/(kgs^2)
    %IERS Table 1.1
    equatorialGravity=9.7803278;%m/s^2
    
    %The astronomical unit. This is roughly the distance from the Earth to
    %the sun.
    AstronomialUnit=1.49597870700*10^11;%meters
    
    %SI standards
    absoluteZero=-273.15;%Centigrade a.k.a. Celsius
    %This is also an ICAO's International Standard Atmosphere's pressure
    standardAtmosphericPressure=101325;%Pascals
    
    %This is the ICAO's International Standard Atmosphere's temperature
    standardTemp=273.15+15;%Degrees Kelvin (15 degrees C)
    %This is the ICAO's International Standard Atmosphere's fraction
    %relative humidity (dry air)
    standardRelHumid=0;
    
    %Time standards
    %GPS time minus international atomic time provides this constant offset
    %in seconds. UTC should not be used in a tracker due to leap seconds
    %interrupting things.
    GPS2TAIOffset=-19;%seconds
    
    %Misc Constants
    %The Euler-Mascheroni constant from mathematics to 128 places. This is
    %the limiting difference between the harmonic series and the natural
    %logarithm. That is
    %EulerMascheroni=lim_{n->Inf}(\sum_{k=1}^n(1/k)-log(n))
    EulerMascheroni=0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951463144724980708248096050401;
    %The golden ratio. Two quantities a and b are in a golden ratio is
    %a/b=(a+b)/a. This is that ratio.
    GoldenRatio=(sqrt(5)+1)/2
    %The Golden Angle. This is the lesser of the two angles that arise from
    %sectioning a circle such that the ratio of the length of the smaller
    %arc to the length of the greater arc is the golden ratio.
    GoldenAngle=pi*(3-sqrt(5));
    %Khinchin's constant to 128 places. This the the geometric mean of the
    %coefficients of almost all continued fractions.
    Khinchin=2.6854520010653064453097148354817956938203822939944629530511523455572188595371520028011411749318476979951534659052880900828976777;
end

methods(Static)
    function [weight,symbol,name]=elementAMU(atomicNumber)
    %%ELEMENTAMU  Get the weight of an element in atomic mass units given
    %             the atmoic number of the element and, if desired, get the
    %             symbol and name of the element as well. If an element has
    %             no stable isotopes, then NaN is returned.
    %
    %INPUTS: atomicNumber An NX1 vector of atomic number of elements of
    %                     interest.
    %
    %OUTPUTS: weight An NX1 vector of the standard atomic weights of the
    %                selected elements in atomic mass units.
    %         symbol An NX1 cell array of text symbols of the elements,
    %                such as 'Na' or 'Au'.
    %         name   An NX1 cell array of text names of the elements.
    %
    %The atomic weights are from the 2013 Table of Standard Atomic Weights
    %that is published by the Commission on Isotopic Abundances and Atomic
    %Weights, Commission II.I of the International Union of Pure and
    %Applied Chemistry and is available at
    %http://www.ciaaw.org/publications.htm
    %Sometimes, the CIAAW gives a range of values for an element instead of
    %just one value. In such an instance, the weight provided is the
    %average of the endpoints of the range.
    %
    %March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.

        dataTable={...
        1, 'hydrogen',        'H',     mean([1.00784, 1.00811]);
        2, 'helium',          'He',    4.002602;
        3, 'lithium',         'Li',    mean([6.938, 6.997]);
        4, 'beryllium',       'Be',    9.0121831;
        5, 'boron',           'B',     mean([10.806, 10.821]);
        6, 'carbon',          'C',     mean([12.0096, 12.0116]);
        7, 'nitrogen',        'N',     mean([14.00643, 14.00728]);
        8, 'oxygen',          'O',     mean([15.99903, 15.99977]);
        9, 'fluorine',        'F',     18.998403163;
        10,'neon',            'Ne',    20.1797;
        11,'sodium',          'Na',    22.98976928;
        12,'magnesium',       'Mg',    mean([24.304, 24.307]);
        13,'aluminium',       'Al',    26.9815385;
        14,'silicon',         'Si',    mean([28.084, 28.086]);
        15,'phosphorus',      'P',     30.973761998;
        16,'sulfur',          'S',     mean([32.059, 32.076]);
        17,'chlorine',        'Cl',    mean([35.446, 35.457]);
        18,'argon',           'Ar',    39.948;
        19,'potassium',       'K',     39.0983;
        20,'calcium',         'Ca',    40.078;
        21,'scandium',        'Sc',    44.955908;
        22,'titanium',        'Ti',    47.867;
        23,'vanadium',        'V',     50.9415;
        24,'chromium',        'Cr',    51.9961;
        25,'manganese',       'Mn',    54.938044;
        26,'iron',            'Fe',    55.845;
        27,'cobalt',          'Co',    58.933194;
        28,'nickel',          'Ni',    58.6934;
        29,'copper',          'Cu',    63.546;
        30,'zinc',            'Zn',    65.38;
        31,'gallium',         'Ga',    69.723;
        32,'germanium',       'Ge',    72.630;
        33,'arsenic',         'As',    74.921595;
        34,'selenium',        'Se',    78.971;
        35,'bromine',         'Br',    mean([79.901, 79.907]);
        36,'krypton',         'Kr',    83.798;
        37,'rubidium',        'Rb',    85.4678;
        38,'strontium',       'Sr',    87.62;
        39,'yttrium',         'Y',     88.90584;
        40,'zirconium',       'Zr',    91.224;
        41,'niobium',         'Nb',    92.90637;
        42,'molybdenum',      'Mo',    95.95;
        43,'technetium',      'Tc',    NaN;
        44,'ruthenium',       'Ru',    101.07;
        45,'rhodium',         'Rh',    102.90550;
        46,'palladium',       'Pd',    106.42;
        47,'silver',          'Ag',    107.8682;
        48,'cadmium',         'Cd',    112.414;
        49,'indium',          'In',    114.818;
        50,'tin',             'Sn',    118.710;
        51,'antimony',        'Sb',    121.760;
        52,'tellurium',       'Te',    127.60;
        53,'iodine',          'I',     126.90447;
        54,'xenon',           'Xe',    131.293;
        55,'caesium',         'Cs',    132.90545196;
        56,'barium',          'Ba',    137.327;
        57,'lanthanum',       'La',    138.90547;
        58,'cerium',          'Ce',    140.116;
        59,'praseodymium',    'Pr',    140.90766;
        60,'neodymium',       'Nd',    144.242;
        61,'promethium',      'Pm',    NaN;
        62,'samarium',        'Sm',    150.36;
        63,'europium',        'Eu',    151.964;
        64,'gadolinium',      'Gd',    157.25;
        65,'terbium',         'Tb',    158.92535; 
        66,'dysprosium',      'Dy',    162.500;
        67,'holmium',         'Ho',    164.93033;
        68,'erbium',          'Er',    167.259;
        69,'thulium',         'Tm',    168.93422;
        70,'ytterbium',       'Yb',    173.054;
        71,'lutetium',        'Lu',    174.9668;
        72,'hafnium',         'Hf',    178.49;
        73,'tantalum',        'Ta',    180.94788;
        74,'tungsten',        'W',     183.84;
        75,'rhenium',         'Re',    186.207;
        76,'osmium',          'Os',    190.23;
        77,'iridium',         'Ir',    192.217;
        78,'platinum',        'Pt',    195.084;
        79,'gold',            'Au',    196.966569;
        80,'mercury',         'Hg',    200.592;
        81,'thallium',        'Tl',    mean([204.382, 204.385]);
        82,'lead',            'Pb',    207.2;
        83,'bismuth',         'Bi',    208.98040; 
        84,'polonium',        'Po',    NaN;
        85,'astatine',        'At',    NaN;
        86,'radon',           'Rn',    NaN;
        87,'francium',        'Fr',    NaN;
        88,'radium',          'Ra',    NaN;
        89,'actinium',        'Ac',    NaN;
        90,'thorium',         'Th',    232.0377;
        91,'protactinium',    'Pa',    231.03588;
        92,'uranium',         'U',     238.02891;
        93,'neptunium',       'Np',    NaN;
        94,'plutonium',       'Pu',    NaN;
        95,'americium',       'Am',    NaN;
        96,'curium',          'Cm',    NaN;
        97,'berkelium',       'Bk',    NaN;
        98,'californium',     'Cf',    NaN;
        99,'einsteinium',     'Es',    NaN;
        100,'fermium',        'Fm',    NaN;
        101,'mendelevium',    'Md',    NaN;
        102,'nobelium',       'No',    NaN;
        103,'lawrencium',     'Lr',    NaN;
        104,'rutherfordium',  'Rf',    NaN;
        105,'dubnium',        'Db',    NaN;
        106,'seaborgium',     'Sg',    NaN;
        107,'bohrium',        'Bh',    NaN;
        108,'hassium',        'Hs',    NaN;
        109,'meitnerium',     'Mt',    NaN;
        110,'darmstadtium',   'Ds',    NaN;
        111,'roentgenium',    'Rg',    NaN;
        112,'copernicium',    'Cn',    NaN;
        113,'ununtrium',      'Uut',   NaN;
        114,'flerovium',      'Fl',    NaN;
        115,'ununpentium',    'Uup',   NaN;
        116,'livermorium',    'Lv',    NaN;
        117,'ununseptium',    'Uus',   NaN;
        118,'ununoctium',     'Uuo',   NaN};

        weight=dataTable(atomicNumber,4);
        weight=cell2mat(weight);
        
        symbol=dataTable(atomicNumber,3);
        name=dataTable(atomicNumber,2);
    end

    function [AMU,C0p,B,dBdT,d2BdT2]=gasProp(substance,T)
    %%GASPROP  Get the atomic weight and, if desired, the ideal gas
    %          specific heat at a constant pressure, and the second virial
    %          coefficient and its is derivatives of various gaseous
    %          substances. Virial coefficients play a role in the physical
    %          chemistry of gasses including the computation of the speed
    %          of sound in a gaseous mixture.
    %
    %INPUTS: substance A text string of the chemical formula of the
    %                  substance. This can be
    %                  'N2'  Nitrogen
    %                  'O2'  Oxygen
    %                  'Ar'  Argon
    %                  'CO2' Carbon dioxide
    %                  'Ne'  Neon
    %                  'Kr'  Krypton
    %                  'CH4' Methane
    %                  'He'  Helium
    %                  'N2O' Nitrous Oxide
    %                  'NO'  Nitrogen Oxide
    %                  'Xe'  Xenon
    %                  'CO'  Carbon Monoxide
    %                  'H2'  Hydrogen
    %                  'H2O' Water
    %                T An NX1 vector of temperatures in degrees Kelvin at
    %                  which the ideal gas specific heat and the second
    %                  virial coefficient of the selected
    %                  substance are  desired. This parameter is not needed
    %                  if only the atomic weight of the substance is
    %                  desired. A warning will be issued if interpolation
    %                  outside of the validated parameter region for the
    %                  substance is performed.
    %
    %OUTPUTS: If an unkown substance is provided, then all empty matrices
    %         will be returned. otherwise:
    %        AMU   The atomic weight of the substance in atomic mass units
    %              (AMU), which is the weight in GRAMS per mol.
    %         C0p  An NX1 vector of the ideal gas specific heat at a
    %              constant pressure for the specific substance at all of
    %              the given temperatures in units of Joules per
    %              kilogram-Kelvin --SI units.
    %         B    An NX1 vector of the second virial coefficients of the
    %              selected substance in SI units of cubic meters per mol
    %              at the provided Kelvin temperatures.
    %       dBdT   An NX1 vector of derivatives of the second virial
    %              coefficient with respect to temperature evaluated at the
    %              selected temperatures, in units of cubic peters per
    %              mol-Kelvin.
    %       d2BdT2 An NX1 vector of second derivatives of the second virial
    %              coefficient with respect to temperature evaluated at the
    %              selected temperatures, in units of cubic peters per
    %              mol-Kelvin^2.
    %
    %The data for all of the second virial coefficients except Helium and
    %water were was taken from
    %D. Ambrose, M. B. Ewing, and M. L. McGlashan. (2014, Mar.) Kaye & Laby
    %tables of physical & chemical constants. National Physical Laboratory.
    %[Online]. Available:
    %http://www.kayelaby.npl.co.uk/chemistry/3_5/3_5.html
    %The second virial coefficient for Helium is taken from cubic spline
    %interpolation of the tabulated values in
    %J. M. H. Levelt Sengers, M. Klein, and J. S.Gallagher, "Pressure-
    %Volume-Temperature Relationships of Gases; Virial Coefficients," in
    %American Institute of Physics Handbook, edited by D. E. Grey
    %(McGraw-Hill, New York, 1972), 3rd ed., Chap. 4i, pp. 4-204-4-221.
    %The second virial coefficient for water was taken from
    %R. W. Hyland, "A correlation for the second interaction virial
    %coefficients and enhancement factors for moist air," Journal of
    %Research of National Bureau of Standards - A Physics and Chemistry,
    %vol. 79A, no. 4, pp. 551-560, Jul. - Aug. 1975.
    %where it was assumed that the valid range matched the temperature
    %range studied in the paper.
    %
    %The data for all of the ideal gas specific heats at a constant
    %pressure are taken from
    %Y. S. Touloukian and T. Makita, "Specific heat nonmetallic liquids and
    %gasses," in Thermophysical Properties of Matter. New York: IFI/Plenum,
    %1970, vol. 6.
    %The tabulated values are in units of cal/(g*K) and thus are converted
    %to SI units.
    %
    %March 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
        abcVirial2=[];
         
        switch(substance)
            case 'N2'%Nitrogen
                AMU=2*Constants.elementAMU(7);
                
                if(nargout>1)
                    C0pRange=[250;1500];

                    C0p=zeros(size(T));
                    sel=T<775;
                    %For T=250K-775K
                    C0p(sel)=0.259934-8.42119e-5*T(sel)+1.72117e-7*T(sel).^2-6.72914e-11*T(sel).^3;
                    %For T= 775K to 1500K
                    C0p(~sel)=0.201678+1.08013e-4*T(~sel)-3.32212e-8*T(~sel).^2+2.45228e-12*T(~sel).^3;
                end
                
                if(nargout>2)
                    abcVirial2=[185.4;141.8;88.7];
                    virialRange=[75;700];
                end
            case 'O2'%Oxygen
                AMU=2*Constants.elementAMU(8);
                
                if(nargout>1)
                    C0pRange=[250;1500];

                    C0p=zeros(size(T));
                    sel=T<760;
                    %For 250K to 760K
                    C0p(sel)=0.222081-7.69230e-5*T(sel)+2.78765e-7*T(sel).^2-1.70107e-10*T(sel).^3;
                    %For 760K to 1500K
                    C0p(~sel)=0.177100+1.49509e-4*T(~sel)-8.44940e-8*T(~sel).^2+1.83236e-11*T(~sel).^3;
                end
                
                if(nargout>2)
                    abcVirial2=[152.8;117.0;108.8];
                    virialRange=[90;400];
                end
            case 'Ar'%Argon
                AMU=Constants.elementAMU(18);
                
                if(nargout>1)
                	C0pRange=[10;6000];
                
                    C0p=ones(size(T))*0.12436;
                end
                
                if(nargout>2)
                    abcVirial2=[154.2;119.3;105.1];
                    virialRange=[80;1024];
                end
            case 'CO2'%Carbon Dioxide
                AMU=Constants.elementAMU(6)+2*Constants.elementAMU(8);
                
                if(nargout>1)
                	C0pRange=[200;1500];

                    C0p=zeros(size(T));
                    sel=T<590;
                    %For T=200K to 590K
                    C0p(sel)=0.105914+4.03552e-4*T(sel)-3.03235e-7*T(sel).^2+8.29431e-11*T(sel).^3;
                    %For T=590K to 1500K
                    C0p(~sel)=0.135069+2.89483e-4*T(~sel)-1.64998e-7*T(~sel)^2+3.53157e-11*T(~sel)^3;
                end
                
                if(nargout>2)
                    abcVirial2=[137.6;87.7;325.7];
                    virialRange=[220;1100];
                end
            case 'Ne'%Neon
                AMU=Constants.elementAMU(10);
                
                if(nargout>1)
                    C0pRange=[10;8000];

                    C0p=ones(size(T))*0.24615;
                end
                
                if(nargout>2)
                    virialRange=[44;973];
                    abcVirial2=[81.0;63.6;30.7];
                end
            case 'Kr'%Krypton
                AMU=Constants.elementAMU(36);
                
                if(nargout>1)
                	C0pRange=[10;6200];
                
                    C0p=ones(size(T))*0.059284;
                end
                
                if(nargout>2)
                    abcVirial2=[189.6;148.0;145.3];
                    virialRange=[110;700];
                end
            case 'CH4'%Methane
                AMU=Constants.elementAMU(6)+4*Constants.elementAMU(1);
                
                if(nargout>1)
                	C0pRange=[270;1500];
                
                    C0p=zeros(size(T));
                    sel=T<790;
                    %For T=270K to 790
                    C0p(sel)=0.458066-2.61341e-4*T(sel)+2.07904e-6*T(sel).^2-1.25017e-9*T(sel).^3;
                    %For T=790K to 1500K
                    C0p(~sel)=0.0258866+1.60802e-3*T(~sel)-6.67069e-7*T(~sel).^2+1.06432e-10*T(~sel).^3;
                end
                
                if(nargout>2)
                    abcVirial2=[206.4;159.5;133.0];
                    virialRange=[110;600];
                end
            case 'He'%Helium
                AMU=Constants.elementAMU(2);

                if(nargout>1)
                    C0pRange=[10;6000];
                    C0p=ones(size(T))*1.2412;
                end
                
                if(nargout>2)
                    %Page 4-206
                    virialTable=[9,     -26.0;
                                 10,    -21.7;
                                 11,    -18.1;
                                 12,    -15.2;
                                 13,    -12.7;
                                 14,    -10.5;
                                 15,    -8.7;
                                 16,    -7.1;
                                 17,    -5.6;
                                 18,    -4.3;
                                 19,    -3.2;
                                 20,    -2.2;
                                 22,    -0.5;
                                 22.64,  0.0;
                                 24,     0.9;
                                 26,     2.0;
                                 28,     3.0;
                                 30,     3.8;
                                 35,     5.4;
                                 40,     6.6;
                                 45,     7.5;
                                 50,     8.2;
                                 60,     9.2;
                                 80,     10.6;
                                 100,    11.4;
                                 120,    11.8;
                                 160,    12.3;
                                 200,    12.3;
                                 273.15, 12.0;
                                 373.15, 11.3;
                                 400,    11.1;
                                 600,    10.4;
                                 800,    9.8;
                                 1000,   9.3;
                                 1200,   8.8;
                                 1400,   8.4];
                    virialRange=[virialTable(1,1);virialTable(end,1)];    
                    
                    %Perform spline interpolation to get the values and
                    %its first two derivatives.
                    BInterp=cubSplineInterpSimp(virialTable(:,1),virialTable(:,2),T);
                    
                    %Convert from cubic centimeters to cubic meters
                    B=BInterp(1,:)'*1e-6;
                    dBdT=BInterp(2,:)'*1e-6;
                    d2BdT2=BInterp(3,:)'*1e-6;
                end
            case 'N2O'%Nitrous Oxide
                AMU=2*Constants.elementAMU(7)+Constants.elementAMU(8);
                
                if(nargout>1)
                    C0pRange=[200;1500];
                
                    C0p=zeros(size(T));
                    sel=T<600;
                    %For T=200K to 600K
                    C0p(sel)=0.103451+4.89293e-4*T(sel)-5.19278e-7*T(sel).^2+2.44839e-10*T(sel).^3;
                    %For T=600K to 1500K
                    C0p(~sel)=0.149343+2.67192e-4*T(~sel)-1.47619e-7*T(~sel).^2+3.01604e-11*T(~sel).^3;
                end
                
                if(nargout>2)
                    abcVirial2=[180.7;114.8;305.4];
                    virialRange=[200;423];
                end
            case 'NO'%Nitrogen Oxide
                AMU=Constants.elementAMU(7)+Constants.elementAMU(8);
                
                if(nargout>1)
                	C0pRange=[100;1500];
                
                    C0p=zeros(size(T));
                    sel=T<590;
                    %For T=100K to 590K
                    C0p(sel)=0.282183-3.16841e-4*T(sel)+6.88734e-7*T(sel).^2-4.22833e-10*T(sel).^3;
                    %For T=590K to 1500K
                    C0p(~sel)=0.192074+1.22287e-4*T(~sel)-5.05602e-8*T(~sel).^2+6.90346e-12*T(~sel).^3;
                end
                
                if(nargout>2)
                    abcVirial2=[15.9;11.0;372.3];
                    virialRange=[122;311];
                end
            case 'Xe'%Xenon
                AMU=Constants.elementAMU(54);
                
                if(nargout>1)
                    C0pRange=[10;5200];

                    C0p=ones(size(T))*0.037837;
                end
                
                if(nargout>2)
                    virialRange=[160;650];
                    abcVirial2=[245.6;190.9;200.2];
                end
            case 'CO'%Carbon Monoxide
                AMU=Constants.elementAMU(6)+Constants.elementAMU(8);
                
                if(nargout>1)
                	C0pRange=[250;1500];
                
                    C0p=zeros(size(T));
                    sel=T<615;
                    %For T=250K to 615K
                    C0p(sel)=0.256859-6.46329e-5*T(sel)+1.31865e-7*T(sel).^2-2.65440e-11*T(sel).^3;
                    %For T=615K to 1500K
                    C0p(~sel)=0.210345+9.44224e-5*T(~sel)-1.94071e-8*T(~sel).^2-2.35385e-12*T(~sel).^3;
                end
                
                if(nargout>2)
                    abcVirial2=[202.6;154.2;94.2];
                    virialRange=[90;573];
                end
            case 'H2'%Hydrogen
                AMU=2*Constants.elementAMU(1);
                
                if(nargout>1)
                	C0pRange=[100;1500];
                
                    C0p=zeros(size(T));
                    sel=T<400;
                    %For T=100K to 400K
                    C0p(sel)=1.46910+1.60057e-2*T(sel)-4.44048e-5*T(sel).^2+4.21220e-8*T(sel).^3;
                    %For T=400K to 1500K
                    C0p(~sel)=3.56903-4.89590e-4*T(~sel)+6.22549e-7*T(~sel).^2-1.19686e-10*T(~sel).^3;
                end
                
                if(nargout>2)
                    abcVirial2=[315.0;289.7;9.47];
                    virialRange=[14;400];
                end
            case 'H2O'%Water
                AMU=2*Constants.elementAMU(1)+Constants.elementAMU(8);
                
                if(nargout>1)
                	C0pRange=[270;1500];
                
                    C0p=zeros(size(T));
                    sel=T<800;
                    %For T=270K to 800K
                    C0p(sel)=0.452219-1.29224e-4*T(sel)+4.17008e-7*T(sel).^2-2.00401e-10*T(sel).^3;
                    %For T=800K to 1500K
                    C0p(~sel)=0.378278+1.53443e-4*T(~sel)+3.31531e-8*T(~sel).^2-1.78435e-11*T(~sel).^3;
                end
                
                if(nargout>2)
                    virialRange=[223.15; 363.15];

                    B=33.97-(55306./T).*10.^(72000./T.^2);
                    dBdT=(27653*2.^(1+72000./T.^2).*5.^(72000./T.^2).*(T.^2+144000*log(10)))./(T.^4);
                    d2BdT2=-((27653*10.^(72000./T.^2).*(T.^4+360000*T.^2*log(10)+10368000000*log(10)^2)*(cosh(2*log(2))+sinh(2*log(2))))./T.^7);
                
                    %Convert from cubic centimeters to cubic meters
                    B=B*1e-6;
                    dBdT=dBdT*1e-6;
                    d2BdT2=d2BdT2*1e-6;
                end
            otherwise%An unknown substance is provided.
                AMU=[];
                C0p=[];
                B=[];
                dBdT=[];
                d2BdT2=[];
                return;
        end
        
        if(nargout>1)
            if(sum(T<C0pRange(1)|T>C0pRange(2))>0)
                warning('A temperature outside of the modelled range (%fK-%fK) for the ideal gas specific heat constant of %s was provided. The results might be unreliable.',C0pRange(1),C0pRange(2),substance);
            end
            
            %Convert from units of Calories per gram-Kelvin to Joules
            %per kilogram-Kelvin.
            %(Calories-> Joules)*(1/g->1/kg)
            C0p=C0p*4.184*1000;
        end
        
        if(nargout>2)
            if(sum(T<virialRange(1)|T>virialRange(2))>0)
                warning('A temperature outside of the modelled range (%fK-%fK) for the second virial coefficient of %s was provided. The results might be unreliable.',virialRange(1),virialRange(2),substance);
            end
        end
        
        %If the second virial coefficient is to be found using a standard
        %exponential model.
        if(~isempty(abcVirial2))
            a=abcVirial2(1);
            b=abcVirial2(2);
            c=abcVirial2(3);

            B=a-b*exp(c./T);
            dBdT=b*c*exp(c./T)./T.^2;
            d2BdT2=-b*c*(c+2*T).*exp(c./T)./T.^4;
            
            %Convert from cubic centimeters to cubic meters
            B=B*1e-6;
            dBdT=dBdT*1e-6;
            d2BdT2=d2BdT2*1e-6;
        end
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
