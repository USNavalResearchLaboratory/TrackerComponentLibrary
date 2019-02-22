function val=convertEnergyUnits(val,unitOrig,unitDest)
%%CONVERTENERGYUNITS Convert values of energy from one set of units to
%                    another. Note that multiple definitions for BTUs cand
%                    calories exist.
%
%INPUTS: val The matrix or vector of values that are to be converted.
%   unitOrig, unitDes Two character strings indicating the units of val and
%            the units into which it is to be converted. The possible
%            values are:
%            'cal'     thermochemical calorie according to ISO
%                      31-4 
%            'calIS'   calorie according to the International
%                      Steam Table of 1956
%            'Cal'     a kilocalorie using the ISO 31-4 definition
%            'CalIS'   a kilocalorie using the International Steam
%                      Table definition
%            'kcal'    synonym for Cal
%            'kcalIS'  synonym for CalIS
%            'BTU'     British Thermal Unit according to ISO 31-4 
%            'BTUPar'  British Thermal Unit according to the
%                      British parliament
%            'BTUIS'   British thermal unit according to the
%                      International Steam Table of 1956 (In the
%                      NIST references).
%            'eV'      electron Volt
%            'erg'     erg (g*cm^2/s^2)
%            'tTNT'    ton of TNT (NIST definition)
%            'ktTNT'   kiloton of TNT
%            'MtTNT'   Megaton of TNT
%            'therm'   Therm as defined for the United States
%            'thermEU' Therm as defined in the European Union
%            'kWh'  kiloWatt-hour 'Wh' Watt-hour
%            'Ws'   Watt-second, synonym for Joule.
%            'YJ'   yottaJoules  'ZJ'  zettaJoules 
%            'EJ'   exaJoules    'PJ'  petaJoules 
%            'TJ'   tetraJoules  'GJ'  gigaJoules 
%            'MJ    megaJoules   'kJ'  kiloJoules 
%            'hJ'   hectoJoules  'daJ' decaJoules 
%            'J'    Joules       'dJ'  deciJoules
%            'cJ'   centiJoules  'mJ'  milliJoules
%            'muJ'  microJoules  'nJ'  nanoJoules
%            'pJ'   picoJoules   'fJ'  femtoJoules
%            'aJ'   attoJoules   'zJ'  zeptoJoules
%            'yJ    yoctoJoules
%
%OUTPUTS: val  The values converted into the desired coordinate system.
%
%%For simplicity, all units are first converted to Joules, and then to the
%desired set of units.
%
%The definitions for many of the conversions are taken from
%http://www.nist.gov/pml/wmd/metric/upload/SP1038.pdf
%http://physics.nist.gov/Pubs/SP811/appenB.html
%http://physics.nist.gov/Pubs/SP811/appenB9.html
%For metric conversions [1] is used.
%
%The definition of the British thermal unit (BTU) from the British
%Parliament's resultation [2]. and lists 1 BTU= 1.05505585257348
%kilojoules. The historical definition of the BTU is according to the
%Oxford English Dictionary online
%( http://www.oed.com ), consulted on 5 May 2014,
%"the amount of heat needed to raise the temperature of one pound of water
%by one degree Fahrenheit..."
%That is a rather vague definition, and is one reason why the BTU should
%not be used as a unit of energy anymore. Similar issues exist with the
%calorie. The conversions used for the Calorie are from the NIST's tables.
%
%The value for the electron volt is taken from the Constants class.
%
%REFERENCES:
%[1] Le Système international d'unités The International System of Units,
%    Bureau international des points et mesures Std., 2006.
%    [Online]. Available: http://www.bipm.org/en/si/si brochure/
%[2] "The Units of Measurement Regulations 1995, 1995 No. 1804" 
%    [Online]. Available: http://www.legislation.gov.uk/uksi/1995/1804/made
%
%May 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%First, convert the value to Joules.
switch(unitOrig)
    case 'cal'%The ISO calorie.
        val=val*4.184;
    case 'calIS'%The calorie from the Fifth International Conference on the
                %Properties ofSteam
        val=val*4.1868;
    case 'Cal'%ISO
        val=val*4.184e3;
    case 'CalIS'%Steam conference
        val=val*4.1868e3;
    case 'kcal'%%Same as Cal
        val=val*4.184e3;
    case 'kcalIS'%Same as CalIS
        val=val*4.1868e3;
    case 'BTU'%BTU as defined by ISO 31-4.
        val=val*1055.056;
    case 'BTUPar'%BTU as defined by the British parliament.
        val=val*1.05505585257348e3;
    case 'BTUIS'%BTU as defined by the Fifth International Conference on
                %the Properties of Steam
        val=val*1.05505585262e3;
    case 'eV'
        val=val*Constants.electronVolt;
    case 'erg'
        val=val*1e-7;
    case 'tTNT'
        val=val*4.184e9;
    case 'ktTNT'
        val=val*4.184e12;
    case 'MtTNT'
        val=val*4.184e15;
    case 'therm'
        val=val*1.054804e8;
    case 'thermEU'
        val=val*1.05506e8;
    case 'kWh'
        val=val*3.6e6;
    case 'Wh'
        val=val*3.6e3;
    case 'Ws'%Same as the Joule
    case 'YJ'
        val=val*1e24;
    case 'ZJ'
        val=val*1e21;
    case 'EJ'
        val=val*1e18;
    case 'PJ'
        val=val*1e15;
    case 'TJ'
        val=val*1e12;
    case 'GJ'
        val=val*1e9;
    case 'MJ'
        val=val*1e6;
    case 'kJ'
        val=val*1e3;
    case 'hJ'
        val=val*1e2;
    case 'daJ'
        val=val*1e1;
    case 'J'%Same as Watt-second
    case 'dJ'
        val=val*1e-1;
    case 'cJ'
        val=val*1e-2;
    case 'mJ'
        val=val*1e-3;
    case 'muJ'
        val=val*1e-6;
    case 'nJ'
        val=val*1e-9;
    case 'pJ'
        val=val*1e-12;
    case 'fJ'
        val=val*1e-15;
    case 'aJ'
        val=val*1e-18;
    case 'zJ'
        val=val*1e-21;
    case 'yJ'
        val=val*1e-24;
    otherwise
        error('The units provided for the source value are invalid.')
end

%Now, convert the value in Joules to the desired coordinate system. The
%multiplications above thus become divisions.
switch(unitDest)
    case 'cal'%The ISO calorie.
        val=val/4.184;
    case 'calIS'%The calorie from the Fifth International Conference on the
                %Properties ofSteam
        val=val/4.1868;
    case 'Cal'%ISO
        val=val/4.184e3;
    case 'CalIS'%Steam conference
        val=val/4.1868e3;
    case 'kcal'%%Same as Cal
        val=val/4.184e3;
    case 'kcalIS'%Same as CalIS
        val=val/4.1868e3;
    case 'BTU'%BTU as defined by ISO 31-4.
        val=val/1055.056;
    case 'BTUPar'%BTU as defined by the British parliament.
        val=val/1.05505585257348e3;
    case 'BTUIS'%BTU as defined by the Fifth International Conference on
                %the Properties of Steam
        val=val/1.05505585262e3;
    case 'eV'
        val=val/Constants.electronVolt;
    case 'erg'
        val=val/1e-7;
    case 'tTNT'
        val=val/4.184e9;
    case 'ktTNT'
        val=val/4.184e12;
    case 'MtTNT'
        val=val/4.184e15;
    case 'therm'
        val=val/1.054804e8;
    case 'thermEU'
        val=val/1.05506e8;
    case 'kWh'
        val=val/3.6e6;
    case 'Wh'
        val=val/3.6e3;
    case 'Ws'%Same as the Joule
    case 'YJ'
        val=val/1e24;
    case 'ZJ'
        val=val/1e21;
    case 'EJ'
        val=val/1e18;
    case 'PJ'
        val=val/1e15;
    case 'TJ'
        val=val/1e12;
    case 'GJ'
        val=val/1e9;
    case 'MJ'
        val=val/1e6;
    case 'kJ'
        val=val/1e3;
    case 'hJ'
        val=val/1e2;
    case 'daJ'
        val=val/1e1;
    case 'J'%Same as Watt-second
    case 'dJ'
        val=val/1e-1;
    case 'cJ'
        val=val/1e-2;
    case 'mJ'
        val=val/1e-3;
    case 'muJ'
        val=val/1e-6;
    case 'nJ'
        val=val/1e-9;
    case 'pJ'
        val=val/1e-12;
    case 'fJ'
        val=val/1e-15;
    case 'aJ'
        val=val/1e-18;
    case 'zJ'
        val=val/1e-21;
    case 'yJ'
        val=val/1e-24;
    otherwise
        error('The units provided for the destination value are invalid.')
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
