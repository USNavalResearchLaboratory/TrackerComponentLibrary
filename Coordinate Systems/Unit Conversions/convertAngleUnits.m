function val=convertAngleUnits(val,unitOrig,unitDest)
%%CONVERTANGLEUNITS Convert values of plane angles from one set of units to
%                   another.
%
%INPUTS: val           The matrix or vector of values that are to be
%                      converted.
%   unitOrig, unitDes  Two character strings indicating the units of val
%                      and the units into which it is to be converted. The
%                      possible values are:
%                      'rad'  radians               'deg' degrees
%                      'amin' arcminute             'as'  arcsecond
%                      'mas'  milliarcsecond        'mus' microarcsecond
%                      'mil'  mil (1/6400 a circle, the NATO usage)
%                      'grad' grad (1/400 a circle) aka gon.
%                      'rev'  revolution (unit of 360 degrees) aka turn.
%                      'quad' quadrant (1/4 a circle)
%                      'sext' sextant (1/6 a circle)
%                      'ha'   hour angle (1/24 a circle)
%                      'pt'   point (1/32 a circle)
%                      'sign' sign (1/12 a circle)
%
%OUTPUTS: val  The values converted into the desired coordinate system.
%
%For simplicity, all units are first converted to radians, and then to the
%desired set of units.
%
%The conversions whose factors are fairely well known and can be found in
%many reference books.
%
%April 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%First, convert the value to radians
switch(unitOrig)
    case 'rad'
    case 'deg'
        val=val*(pi/180);
    case 'amin'
        %arcminutes->degrees->radians
        val=val*((1/60)*(pi/180));
    case 'as'
        %arcseconds->arcminutes->degrees->radians
        val=val*((1/60)*(1/60)*(pi/180));
    case 'mas'
        %milliarcseconds->arcseconds->arcminutes->degrees->radians
        val=val*(1e-3*(1/60)*(1/60)*(pi/180));
    case 'mus'
        %microarcseconds->arcseconds->arcminutes->degrees->radians
        val=val*(1e-6*(1/60)*(1/60)*(pi/180));
    case 'mil'
        val=val*(2*pi/6400);
    case 'grad'
        val=val*(2*pi/400);
    case 'rev'
        val=val*(2*pi);
    case 'quad'
        val=val*(2*pi/4);
    case 'sext'
        val=val*(2*pi/6);
    case 'ha'
        val=val*(2*pi/24);
    case 'pt'
        val=val*(2*pi/32);
    case 'sign'
        val=val*(2*pi/12);
    otherwise
        error('The units provided for the source value are invalid.')
end

%Now, convert the value in radians to the desired coordinate system. The
%multiplications above thus become divisions.
switch(unitDest)
    case 'rad'
    case 'deg'
        val=val/(pi/180);
    case 'amin'
        %arcminutes->degrees->radians
        val=val/((1/60)*(pi/180));
    case 'as'
        %arcseconds->arcminutes->degrees->radians
        val=val/((1/60)*(1/60)*(pi/180));
    case 'mas'
        %milliarcseconds->arcseconds->arcminutes->degrees->radians
        val=val/(1e-3*(1/60)*(1/60)*(pi/180));
    case 'mus'
        %microarcseconds->arcseconds->arcminutes->degrees->radians
        val=val/(1e-6*(1/60)*(1/60)*(pi/180));
    case 'mil'
        val=val/(2*pi/6400);
    case 'grad'
        val=val/(2*pi/400);
    case 'rev'
        val=val/(2*pi);
    case 'quad'
        val=val/(2*pi/4);
    case 'sext'
        val=val/(2*pi/6);
    case 'ha'
        val=val/(2*pi/24);
    case 'pt'
        val=val/(2*pi/32);
   case 'sign'
        val=val/(2*pi/12);
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
