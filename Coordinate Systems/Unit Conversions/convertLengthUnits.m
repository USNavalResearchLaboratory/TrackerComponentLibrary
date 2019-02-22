function val=convertLengthUnits(val,unitOrig,unitDest)
%%CONVERTLENGTHUNITS Convert values of length from one set of units to
%                    another.
%
%INPUTS: val The matrix or vector of values that are to be converted.
% unitOrig, unitDes Two character strings indicating the units of val and
%            the units into which it is to be converted. The possible
%            values are:
%            'Ym'   yottameters        'Zm'  zettameters
%            'Em'   exameters          'Pm'  petameters
%            'Tm'   tetrameters        'Gm'  gigameters
%            'Mm    megameters         'km'  kilometers
%            'hm'   hectometers        'dam' decameters
%            'm'    meters             'dm'  decimeters
%            'cm'   centimeters        'mm'  millimeters
%            'mum'  micrometers        'nm'  nanometers
%            'pm'   picometers         'fm'  femtometers
%            'am'   attometers         'zm'  zeptometers
%            'ym'   yoctometers        'nml' nautical miles
%            'mi'   miles              'yd'  yards
%            'ft'   feet               'in'  inches
%            'pt'   points             'pica'picas (1/6 inch)
%            'ftm'  fathoms            'AU'  astronomical units
%            'ly'   light-years        'pc'  parsec
%            'Ang'  Angströms          'mil' mils (1/1000 inches)
%            'cbl'  Cable length (0.1nml) Note that many other definitions
%                   of a cable length exist.                
%            'USft' US survey feet     'USmi' US survey miles
%            'USftm'US survey fathoms  
%
%OUTPUTS: val  The values converted into the desired coordinate system.
%
%For simplicity, all units are first converted to meters, and then to the
%desired set of units. The US survey units stem from the fact that some
%maps and offocial documents in the United States kept the old definition
%of the foot and mile rather than switching to newer, standard definitions.
%The US survey mile is also known as a statutory mile. The definition of
%the foot and mile that are commonly used everyday in the United State are
%not the survey values.
%
%The definitions for the conversions are taken from
%For the parsec:
%http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_2.html
%for the astronomical unit
%http://www.iau.org/static/resolutions/IAU2012_English.pdf
%For misc conversions
%http://www.nist.gov/pml/wmd/metric/upload/SP1038.pdf
%http://physics.nist.gov/Pubs/SP811/appenB.html
%http://physics.nist.gov/Pubs/SP811/appenB9.html
%For metric conversions, [1] was used.
%
%REFERENCES:
%[1] Le Système international d'unités The International System of Units,
%    Bureau international des points et mesures Std., 2006.
%    [Online]. Available: http://www.bipm.org/en/si/si brochure/
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%First, convert the value to meters
switch(unitOrig)
    case 'Ym'
        val=val*1e24;
    case 'Zm'
        val=val*1e21;
    case 'Em'
        val=val*1e18;
    case 'Pm'
        val=val*1e15;
    case 'Tm'
        val=val*1e12;
    case 'Gm'
        val=val*1e9;
    case 'Mm'
        val=val*1e6;
    case 'km'
        val=val*1e3;
    case 'hm'
        val=val*1e2;
    case 'dam'
        val=val*1e1;
    case 'm'
    case 'dm'
        val=val*1e-1;
    case 'cm'
        val=val*1e-2;
    case 'mm'
        val=val*1e-3;
    case 'mum'
        val=val*1e-6;
    case 'nm'
        val=val*1e-9;
    case 'pm'
        val=val*1e-12;
    case 'fm' 
        val=val*1e-15;
    case 'am'
        val=val*1e-18;
    case 'zm'
        val=val*1e-21;
    case 'ym'
        val=val*1e-24;
    case 'nml'
        val=val*1.852e3;%NIST
    case 'mi'
        val=val*1.609344e3;%NIST
    case 'yd'
        val=val*0.9144;%NIST
    case 'ft'
        val=val*0.3048;%NIST
    case 'in'
        val=val*0.0254;%NIST.
    case 'pt'
        %points to inches (via w3) and then inches to meters (via NIST).
        val=val*((1/72)*0.0254);
    case 'pica'
        %pica->points->inches->meters
        val=val*(12*(1/72)*0.0254);
    case 'ftm'
        %fathom->feet(From Oxford American Dictionary)->meters (NIST)
        val=val*(6*0.3048);
    case 'AU'
        val=val*149597870700;%From the IAU.
    case 'ly'
        val=val*9460730472580800;
    case 'pc'
        arcSecInRad=(1/60)*(1/60)*(pi/180);
        %A parsec is defined as the distance at which 1 AU subtends 1
        %second of arc. Thus, using the definition of the AU:
        val=val*((1/sin(arcSecInRad))*149597870700);
    case 'Ang'
        val=val*1e-10;%NIST
    case 'mil'
        val=val*2.54e-5;%NIST
    case 'cbl'
        %NIST for nautical mile and Oxford English dictionary for the cable
        %length in terms of nautical miles.
        val=val*1.852e3*0.1;
    case 'USft'%NIST
        val=val*(1200/3937);
    case 'USmi'%NIST
        %US survey miles (statute miles) to US survey feet to meters.
        val=val*(5280*(1200/3937));
    case 'USftm'%NIST
        %Fathoms in US Survey feet to US survey feet to meters.
        val=val*(6*(1200/3937));
    otherwise
        error('The units provided for the source value are invalid.')
end

%Now, convert the value in meters to the desired coordinate system. The
%multiplications above thus become divisions.
switch(unitDest)
    case 'Ym'
        val=val/1e24;
    case 'Zm'
        val=val/1e21;
    case 'Em'
        val=val/1e18;
    case 'Pm'
        val=val/1e15;
    case 'Tm'
        val=val/1e12;
    case 'Gm'
        val=val/1e9;
    case 'Mm'
        val=val/1e6;
    case 'km'
        val=val/1e3;
    case 'hm'
        val=val/1e2;
    case 'dam'
        val=val/1e1;
    case 'm'
    case 'dm'
        val=val/1e-1;
    case 'cm'
        val=val/1e-2;
    case 'mm'
        val=val/1e-3;
    case 'mum'
        val=val/1e-6;
    case 'nm'
        val=val/1e-9;
    case 'pm'
        val=val/1e-12;
    case 'fm' 
        val=val/1e-15;
    case 'am'
        val=val/1e-18;
    case 'zm'
        val=val/1e-21;
    case 'ym'
        val=val/1e-24;
    case 'nml'
        val=val/1.852e3;
    case 'mi'
        val=val/1.609344e3;
    case 'yd'
        val=val/0.9144;
    case 'ft'
        val=val/0.3048;
    case 'in'
        val=val/0.0254;
    case 'pt'
        val=val/((1/72)*0.0254);
    case 'pica'
        val=val/(12*(1/72)*0.0254);
    case 'ftm'
        val=val/(6*0.3048);
    case 'AU'
        val=val/149597870700;
    case 'ly'
        val=val/9460730472580800;
    case 'pc'
        arcSecInRad=(1/60)*(1/60)*(pi/180);
        val=val/((1/sin(arcSecInRad))*149597870700);
    case 'Ang'
        val=val/1e-10;
    case 'mil'
        val=val/2.54e-5;
    case 'cbl'
        val=val/(1.852e3*0.1);
    case 'USft'
        val=val/(1200/3937);
    case 'USmi'%NIST
        val=val/(5280*(1200/3937));
    case 'USftm'%NIST
        val=val/(6*(1200/3937));
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
