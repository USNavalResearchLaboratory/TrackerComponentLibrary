function val=convertTimeUnits(val,unitOrig,unitDest)
%%CONVERTTIMEUNITS Convert durations of time from one set of units to
%                  another. This is for converting durations of time. To
%                  convert actual dates in various timescales, look at the
%                  functions in ./Coordinate Systems/Time .
%
%INPUTS: val The matrix or vector of values that are to be converted.
% unitOrig, unitDes Two character strings indicating the units of val
%            and the units into which it is to be converted. The possible
%            values are:
%            'wk'   weeks         'd'   days
%            'h'    hour          'min' minutes
%            'Ys'   yottaseconds  'Zs'  zettaseconds 
%            'Es'   exaseconds    'Ps'  petaseconds 
%            'Ts'   tetraseconds  'Gs'  gigaseconds 
%            'Ms    megaseconds   'ks'  kiloseconds 
%            'hs'   hectoseconds  'das' decaseconds 
%            's'    seconds       'ds'  deciseconds
%            'cs'   centiseconds  'ms'  milliseconds
%            'mus'  microseconds  'ns'  nanoseconds
%            'ps'   picoseconds   'fs'  femtoseconds
%            'as'   attoseconds   'zs'  zeptoseconds
%            'ys    yoctoseconds  'fn'  fortnights (14 days)
%
%OUTPUTS: val The values converted into the desired coordinate system.
%
%%For simplicity, all units are first converted to seconds, and then to the
%desired set of units.
%
%The metric conversions are taken from [1]. Other conversions tend to be
%common knowledge.
%
%Many intervals such as years and decades are not included, since the
%duration of a year is not  a constant number of seconds.
%
%REFERENCES:
%[1] Le Système international d'unités The International System of Units,
%    Bureau international des points et mesures Std., 2006.
%    [Online]. Available: http://www.bipm.org/en/si/si brochure/
%
%May 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%First, convert the value to seconds.
switch(unitOrig)
    case 'wk'
        %weeks->days->hours->minuts->seconds
        val=val*(7*24*60*60);
    case 'd'
        %days->hours->minuts->seconds
        val=val*(24*60*60);
    case 'h'
        %hours->minuts->seconds
        val=val*(60*60);
    case 'min'
        val=val*60;
    case 'Ys'
        val=val*1e24;
    case 'Zs'
        val=val*1e21;
    case 'Es'
        val=val*1e18;
    case 'Ps'
        val=val*1e15;
    case 'Ts'
        val=val*1e12;
    case 'Gs'
        val=val*1e9;
    case 'Ms'
        val=val*1e6;
    case 'ks'
        val=val*1e3;
    case 'hs'
        val=val*1e2;
    case 'das'
        val=val*1e1;
    case 's'
    case 'ds'
        val=val*1e-1;
    case 'cs'
        val=val*1e-2;
    case 'ms'
        val=val*1e-3;
    case 'mus'
        val=val*1e-6;
    case 'ns'
        val=val*1e-9;
    case 'ps'
        val=val*1e-12;
    case 'fs'
        val=val*1e-15;
    case 'as'
        val=val*1e-18;
    case 'zs'
        val=val*1e-21;
    case 'ys'
        val=val*1e-24;
    case 'fn'
        %fortnights->days->hours->minuts->seconds
        val=val*(14*24*60*60);
    otherwise
        error('The units provided for the source value are invalid.')
end

%Now, convert the value in seconds to the desired coordinate system. The
%multiplications above thus become divisions.
switch(unitDest)
    case 'wk'
        %weeks->days->hours->minuts->seconds
        val=val/(7*24*60*60);
    case 'd'
        %days->hours->minuts->seconds
        val=val/(24*60*60);
    case 'h'
        %hours->minuts->seconds
        val=val/(60*60);
    case 'min'
        val=val/60;
    case 'Ys'
        val=val/1e24;
    case 'Zs'
        val=val/1e21;
    case 'Es'
        val=val/1e18;
    case 'Ps'
        val=val/1e15;
    case 'Ts'
        val=val/1e12;
    case 'Gs'
        val=val/1e9;
    case 'Ms'
        val=val/1e6;
    case 'ks'
        val=val/1e3;
    case 'hs'
        val=val/1e2;
    case 'das'
        val=val/1e1;
    case 's'
    case 'ds'
        val=val/1e-1;
    case 'cs'
        val=val/1e-2;
    case 'ms'
        val=val/1e-3;
    case 'mus'
        val=val/1e-6;
    case 'ns'
        val=val/1e-9;
    case 'ps'
        val=val/1e-12;
    case 'fs'
        val=val/1e-15;
    case 'as'
        val=val/1e-18;
    case 'zs'
        val=val/1e-21;
    case 'ys'
        val=val/1e-24;
    case 'fn'
        %fortnights->days->hours->minuts->seconds
        val=val/(14*24*60*60);
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
