function val=convertForceUnits(val,unitOrig,unitDes)
%%CONVERTFORCEUNITS Convert values of force from one set of units to
%                   another.
%
%INPUTS: val The matrix or vector of values that are to be converted.
% unitOrig, unitDes Two character strings indicating the units of val and
%            the units into which it is to be converted. The possible
%            values are:
%            'dyn'    Dynes (cg*g/s^2)
%            'kip'    1000 pounds-force
%            'kgf'    kilograms-force
%            'lbf'    pounds-force
%            'ozf'    ounces-force
%            'YPa'    yottaNewtons      'ZPa'  zettaNewtons
%            'EPa'    exaNewtons        'PPa'  petaNewtons
%            'TPa'    tetraNewtons      'GPa'  gigaNewtons
%            'MPa'    megaNewtons       'kPa'  kiloNewtons
%            'hPa'    hectoNewtons      'daPa' decaNewtons
%            'Pa'     Newtons(kg*m/s^2) 'dPa'  deciNewtons
%            'cPa'    centiNewtons      'mPa'  milliNewtons
%            'muPa'   microNewtons      'nPa'  nanoNewtons
%            'pPa'    picoNewtons       'fPa'  femtoNewtons
%            'aPa'    attoNewtons       'zPa'  zeptoNewtons
%            'yPa'    yoctoNewtons
%            'pdl'    poundals
%            'tonf'   ton-force, (2000lbf, Note: other definitions exist)
%
%Many of the conversions were taken from [1].
%
%REFERENCES:
%[1] Guide for the Use of the International System of Units, NIST Special
%    Publication 811, National Institute of standards and Technology, 2008
%    Edition.
%    https://www.wmo.int/pages/prog/gcos/documents/gruanmanuals/NIST/sp811.pdf
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%First, convert the value to Newtons
switch(unitOrig)
    case 'dyn'
        val=val*1e-5;
    case 'kip'%1kip=1000lbf
        val=val*(4.4482216152605*1000);
    case 'kgf'
        val=val*9.80665;
    case 'lbf'%pound-force
        val=val*4.4482216152605;
    case 'ozf'%ounce-force (16 ounces per pound)
        val=val*(4.4482216152605*16);
    case 'pdl'%One poundal is the force that accelerates a 1-pound 
          %(Imperial) mass at 1 foot per second^2
        val=val*0.138254954376;
    case 'tonf'%ton force (1 tonf=2000lbf)
        val=val*(4.4482216152605*2000);
    case 'YN'
        val=val*10^(24);
    case 'ZN'
        val=val*10^(21);
    case 'EN'
        val=val*10^(18);
    case 'PN'
        val=val*10^(15);
    case 'TN'
        val=val*10^(12);
    case 'GN'
        val=val*10^(9);
    case 'MN'
        val=val*10^(6);
    case 'kN'
        val=val*10^(3);
    case 'hN'
        val=val*10^(2);
    case 'daN'
        val=val*10^(1);
    case 'N'
    case 'dN'
        val=val*10^(-1);
    case 'cN'
        val=val*10^(-2);
    case 'mN'
        val=val*10^(-3);
    case 'muN'
        val=val*10^(-6);
    case 'nN'
        val=val*10^(-9);
    case 'pN'
        val=val*10^(-12);
    case 'fN'
        val=val*10^(-15);
    case 'aN'
        val=val*10^(-18);
    case 'zN'
        val=val*10^(-21);
    case 'yN'
        val=val*10^(-24);
    otherwise
        error('The units provided for the source value are invalid.')
end

%Now, convert the value in pascals to the desired coordinate system. The
%multiplications above thus become divisions.
switch(unitDes)
    case 'dyn'
        val=val/1e-5;
    case 'kip'%1kip=1000lbf
        val=val/(4.4482216152605*1000);
    case 'kgf'
        val=val/9.80665;
    case 'lbf'%pound-force
        val=val/4.4482216152605;
    case 'ozf'%ounce-force (16 ounces per pound)
        val=val/(4.4482216152605*16);
    case 'pdl'%One poundal is the force that accelerates a 1-pound 
          %(Imperial) mass at 1 foot per second^2
        val=val/0.138254954376;
    case 'tonf'%ton force (1 tonf=2000lbf)
        val=val/(4.4482216152605*2000);
    case 'YN'
        val=val/10^(24);
    case 'ZN'
        val=val/10^(21);
    case 'EN'
        val=val/10^(18);
    case 'PN'
        val=val/10^(15);
    case 'TN'
        val=val/10^(12);
    case 'GN'
        val=val/10^(9);
    case 'MN'
        val=val/10^(6);
    case 'kN'
        val=val/10^(3);
    case 'hN'
        val=val/10^(2);
    case 'daN'
        val=val/10^(1);
    case 'N'
    case 'dN'
        val=val/10^(-1);
    case 'cN'
        val=val/10^(-2);
    case 'mN'
        val=val/10^(-3);
    case 'muN'
        val=val/10^(-6);
    case 'nN'
        val=val/10^(-9);
    case 'pN'
        val=val/10^(-12);
    case 'fN'
        val=val/10^(-15);
    case 'aN'
        val=val/10^(-18);
    case 'zN'
        val=val/10^(-21);
    case 'yN'
        val=val/10^(-24);
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
