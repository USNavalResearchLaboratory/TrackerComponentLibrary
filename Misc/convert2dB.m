function val=convert2dB(val,units,refVal)
%%CONVERT2DB Given a real or complex voltage or a power value, convert the
%            value to decibels (dB).
%
%INPUTS: val A scalar or matrix of values to convert to decibels. These can
%            be complex.
%      units The scaling of the conversion depends on whether a power or
%            voltage is given. Possible values are:
%            'voltage" (The default if omitted or an empty matrix is
%                       passed.)
%            'power'
%     refVal The real reference value that forms the zero-dB level. If val
%            holds powers, then this will typically have units of Watts. If
%            vals hold voltages, then this will typically have units of
%            Volts. If this parameter is omitted or an empty matrix is
%            passed, the default of 1 is used.
%
%OUTPUTS: val The input value converted to decibels.
%
%Decibels are defined in terms of power. If a voltage is given, then the
%output is 10*log10(abs(val).^2./refVal). Otherwise, the output when given
%a power is 10*log10(abs(val)./refVal).
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(refVal))
    refVal=1;
end

if(nargin<2||isempty(units))
    units='voltage'; 
end

switch(units)
    case 'voltage'
        val=10*log10(abs(val).^2./refVal);
    case 'power'
        val=10*log10(abs(val)./refVal);
    otherwise
        error('Unknown units specified')
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
