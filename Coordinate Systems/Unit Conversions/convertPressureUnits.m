function val=convertPressureUnits(val,unitOrig,unitDes)
%%CONVERTPRESSUREUNITS Convert values of pressure from one set of units to
%                      another.
%
%INPUTS: val The matrix or vector of values that are to be converted.
% unitOrig, unitDes Two character strings indicating the units of val and
%            the units into which it is to be converted. The possible
%            values are:
%            'at'  Technical atmospheric units (1kg force/cm^2)
%            'atm' Standard atmospheric units (101325Pa)
%            'Ybar'   yottabar        'Zbar'  zettabar
%            'Ebar'   exabar          'Pbar'  petabar
%            'Tbar'   tetrabar        'Gbar'  gigabar
%            'Mbar'   megabar         'kbar'  kilobar
%            'hbar'   hectobar        'dabar' decabar
%            'bar'    bar 10^6 Pa     'dbar'  decibar
%            'cbar'   centibar        'mbar'  millibar
%            'mubar'  microbar        'nbar'  nanobar
%            'pbar'   picobar         'fbar'  femtobar
%            'abar'   attobar         'zbar'  zeptobar
%            'ybar'   yoctobar
%            'Ba'     Barye 1 dyne/cm^2
%            'inHg'   inches of mercury
%            'mmHg'   millimeters of mercury (1Torr)
%            'YPa'    yottaPascals    'ZPa'  zettaPascals
%            'EPa'    exaPascals      'PPa'  petaPascals
%            'TPa'    tetraPascals    'GPa'  gigaPascals
%            'MPa'    megaPascals     'kPa'  kiloPascals
%            'hPa'    hectoPascals    'daPa' decaPascals
%            'Pa'     Pascals (N/m^2) 'dPa'  deciPascals
%            'cPa'    centiPascals    'mPa'  milliPascals
%            'muPa'   microPascals    'nPa'  nanoPascals
%            'pPa'    picoPascals     'fPa'  femtoPascals
%            'aPa'    attoPascals     'zPa'  zeptoPascals
%            'yPa'    yoctoPascals 
%            'psi'    Pound (Avoirdupois) per square inch                      
%
%1 atmosphere pressure is defined as 101325Pa in
%http://www.bipm.org/fr/CGPM/db/10/4/
%
%The pressure that raises 1mmHg varies depending on the temperature and
%other factors. Thus, the standardized value legally defined by the EU in
%http://eur-lex.europa.eu/LexUriServ/LexUriServ.do?uri=CONSLEG:1980L0181:20090527:DE:PDF
%is used. There 1mmHg is 133.322Pa.
%
%The value for the conversion from pounds to Newtons to convert PSI to
%Pascals is taken from [1]. That is also where the conversion from kilogram
%-force to Newtons was taken and where the definition of the technical
%atmosphere can be found.
%
%REFERENCES:
%[1] Guide for the Use of the International System of Units, NIST Special
%    Publication 811, National Institute of standards and Technology, 2008
%    Edition.
%    https://www.wmo.int/pages/prog/gcos/documents/gruanmanuals/NIST/sp811.pdf
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%First, convert the value to Pascals
switch(unitOrig)
    case 'at'
        %The first term converts from kilogram-force to Newtons. The second
        %term converts from squared centimeters in the denominator to
        %square meters.
        val=val*((9.80665)/(0.0001));
    case 'atm'
        val=val*101325;
    case 'Ybar'
        val=val*10^(5+24);
    case 'Zbar'
        val=val*10^(5+21);
    case 'Ebar'
        val=val*10^(5+18);
    case 'Pbar'
        val=val*10^(5+15);
    case 'Tbar'
        val=val*10^(5+12);
    case 'Gbar'
        val=val*10^(5+9);
    case 'Mbar'
        val=val*10^(5+6);
    case 'kbar'
        val=val*10^(5+3);
    case 'hbar'
        val=val*10^(5+2);
    case 'dabar'
        val=val*10^(5+1);
    case 'bar'
        val=val*10^5;
    case 'dbar'
        val=val*10^(5-1);
    case 'cbar'
        val=val*10^(5-2);
    case 'mbar'
        val=val*10^(5-3);
    case 'mubar'
        val=val*10^(5-6);
    case 'nbar'
        val=val*10^(5-9);
    case 'pbar'
        val=val*10^(5-12);
    case 'fbar'
        val=val*10^(5-15);
    case 'abar'
        val=val*10^(5-18);
    case 'zbar'
        val=val*10^(5-21);
    case 'ybar'
        val=val*10^(5-24);
    case 'Ba'
        val=val/10;
    case 'inHg'
        %The 25.4 converts from inches to millimeters, then the standard
        %for mmHg is used.
        val=val*(25.4*133.322);
    case 'mmHg'
        val=val*133.322;
    case 'YPa'
        val=val*10^(24);
    case 'ZPa'
        val=val*10^(21);
    case 'EPa'
        val=val*10^(18);
    case 'PPa'
        val=val*10^(15);
    case 'TPa'
        val=val*10^(12);
    case 'GPa'
        val=val*10^(9);
    case 'MPa'
        val=val*10^(6);
    case 'kPa'
        val=val*10^(3);
    case 'hPa'
        val=val*10^(2);
    case 'daPa'
        val=val*10^(1);
    case 'Pa'
    case 'dPa'
        val=val*10^(-1);
    case 'cPa'
        val=val*10^(-2);
    case 'mPa'
        val=val*10^(-3);
    case 'muPa'
        val=val*10^(-6);
    case 'nPa'
        val=val*10^(-9);
    case 'pPa'
        val=val*10^(-12);
    case 'fPa'
        val=val*10^(-15);
    case 'aPa'
        val=val*10^(-18);
    case 'zPa'
        val=val*10^(-21);
    case 'yPa'
        val=val*10^(-24);
    case 'psi'
        %First term converts pounds to Newtons. The second term converts
        %square inches in the denominator to square meters, giving Pascals.
        val=val*((4.4482216152605)/(0.0254)^2);
    otherwise
        error('The units provided for the source value are invalid.')
end

%Now, convert the value in pascals to the desired coordinate system. The
%multiplications above thus become divisions.
switch(unitDes)
    case 'at'
        %The first term converts from kilogram-force to Newtons. The second
        %term converts from squared centimeters in the denominator to
        %square meters.
        val=val/((9.80665)/(0.0001));
    case 'atm'
        val=val/101325;
    case 'Ybar'
        val=val/10^(5+24);
    case 'Zbar'
        val=val/10^(5+21);
    case 'Ebar'
        val=val/10^(5+18);
    case 'Pbar'
        val=val/10^(5+15);
    case 'Tbar'
        val=val/10^(5+12);
    case 'Gbar'
        val=val/10^(5+9);
    case 'Mbar'
        val=val/10^(5+6);
    case 'kbar'
        val=val/10^(5+3);
    case 'hbar'
        val=val/10^(5+2);
    case 'dabar'
        val=val/10^(5+1);
    case 'bar'
        val=val/10^5;
    case 'dbar'
        val=val/10^(5-1);
    case 'cbar'
        val=val/10^(5-2);
    case 'mbar'
        val=val/10^(5-3);
    case 'mubar'
        val=val/10^(5-6);
    case 'nbar'
        val=val/10^(5-9);
    case 'pbar'
        val=val/10^(5-12);
    case 'fbar'
        val=val/10^(5-15);
    case 'abar'
        val=val/10^(5-18);
    case 'zbar'
        val=val/10^(5-21);
    case 'ybar'
        val=val/10^(5-24);
    case 'Ba'
        val=val/10;
    case 'inHg'
        %The 25.4 converts from inches to millimeters, then the standard
        %for mmHg is used.
        val=val/(25.4*133.322);
    case 'mmHg'
        val=val/133.322;
    case 'YPa'
        val=val/10^(24);
    case 'ZPa'
        val=val/10^(21);
    case 'EPa'
        val=val/10^(18);
    case 'PPa'
        val=val/10^(15);
    case 'TPa'
        val=val/10^(12);
    case 'GPa'
        val=val/10^(9);
    case 'MPa'
        val=val/10^(6);
    case 'kPa'
        val=val/10^(3);
    case 'hPa'
        val=val/10^(2);
    case 'daPa'
        val=val/10^(1);
    case 'Pa'
    case 'dPa'
        val=val/10^(-1);
    case 'cPa'
        val=val/10^(-2);
    case 'mPa'
        val=val/10^(-3);
    case 'muPa'
        val=val/10^(-6);
    case 'nPa'
        val=val/10^(-9);
    case 'pPa'
        val=val/10^(-12);
    case 'fPa'
        val=val/10^(-15);
    case 'aPa'
        val=val/10^(-18);
    case 'zPa'
        val=val/10^(-21);
    case 'yPa'
        val=val/10^(-24);
    case 'psi'
        %First term converts pounds to Newtons. The second term converts
        %square inches in the denominator to square meters, giving Pascals.
        val=val/((4.4482216152605)/(0.0254)^2);
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
