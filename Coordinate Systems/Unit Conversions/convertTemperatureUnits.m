function val=convertTemperatureUnits(val,unitOrig,unitDest)
%%CONVERTTEMPERATUREUNITS Convert values of temperature from one set of
%                         units to another.
%
%INPUTS: val           The matrix or vector of values that are to be
%                      converted.
%   unitOrig, unitDes  Two character strings indicating the units of val
%                      and the units into which it is to be converted. The
%                      possible values are:
%                      'K' Kelvin    'C' Centigrade/ Celcius
%                      'R' Rankine   'F' Fahrenheit
%
%OUTPUTS: val  The values converted into the desired coordinate system.
%
%The conversions are standard and can be found in many places books and
%places online, such as at
%http://scienceworld.wolfram.com/physics/Rankine.html
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%First, convert the value to Kelvin
switch(unitOrig)
    case 'K'
    case 'C'
        val=val+273.15;
    case 'R'
        val=val*(5/9);
    case 'F'
        val=(val-32)*(5/9)+273.15;
    otherwise
        error('The units provided for the source value are invalid.')
end

%Then convert the temperature in Kelvin to the desired temperature.
switch(unitDest)
    case 'K'
    case 'C'
        val=val-273.15;
    case 'R'
        val=val/(5/9);
    case 'F'
        val=32+(9/5)*(val-273.15);
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

