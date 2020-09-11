function val=convertFromdB(dBVal,units,refVal)
%%CONVERTFROMDB Given a voltage or power value in decibels, convert the
%               value back to a non-decibel value. This function is the
%               opposite of convert2dB.
%
%INPUTS: dBVal A scalar or matrix of real values in decibels.
%       units The scaling of the conversion depends on whether a power or
%             voltage is given. Possible values are:
%             'voltage" (The default if omitted or an empty matrix is
%                        passed.)
%             'power'
%      refVal The postivie, real reference value that forms the zero-dB
%             level. If dBVal holds powers, then this will typically have
%             units of Watts. If dBVal hold voltages, then this will
%             typically have units of Volts. If this parameter is omitted
%             or an empty matrix is passed, the default of 1 is used.
%
%OUTPUTS: val The positive real value
%
%Decibels are defined in terms of power. If a voltage is given, then the
%output is sqrt(refVal.*10.^(dBVal/10)). Otherwise, the output when given
%a power is refVal.*10.^(dBVal/10).
%
%February 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(refVal))
    refVal=1;
end

if(nargin<2||isempty(units))
    units='voltage'; 
end

switch(units)
    case 'voltage'
        val=sqrt(refVal.*10.^(dBVal/10));
    case 'power'
        val=refVal.*10.^(dBVal/10);
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
