classdef RadarBands
%%RADARBANDS A collection of functions for converting between radar
%            frequencies and band name conventions.
%
%Implemented conventions are: ITU, IEEE, NATO
%
%REFERENCES:
%[1] I. T. U. (ITU), "Nomenclature of the frequency and wavelength bands
%    used in telecommunications," Recommendation ITU/RV, pp. 431-438, 2015.
%[2] J. A. Bruder, "IEEE radar standards and the radar systems panel,"
%    IEEE Aerospace and Electronic Systems Magazine, vol. 28, no. 7, pp.
%    19-22, 2013. doi: 10.1109/MAES.2013.6559377.
%[3] L. A. Belov, S. M. Smolskiy, and V. N. Kochemasov, Handbook of RF,
%    microwave, and millimeter-wave components. Artech house, 2012.
%
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

methods(Static)
    function freq = getITUBounds(band)
        %%GETITUBOUNDS Returns the frequency bounds which define the given
        %              ITU frequency band.
        %
        %INPUT:
        % band: A string or character array which identifies the
        %       abbreviated ITU frequency band.
        %
        %OUTPUT:
        % freq: A 1-by-2 vector containing the lower and upper bound for
        %       the specified frequency band specified in Hertz. If band
        %       does not match one of the known abbreviations, [NaN,NaN]
        %       is returned.
        %
        %December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        %
        bounds = [0,3*10.^(0:12)];
        bands = ["TLF","ELF","SLF","ULF","VLF","LF","MF","HF","VHF","UHF","SHF","EHF","THF"];
        i = find(bands==band);
        if any(i)
            freq = bounds(i:i+1);
        else
            freq = [NaN,NaN];
        end
    end
    function freq = getIEEEBounds(band)
        %%GETIEEEBOUNDS Returns the frequency bounds which define the given
        %               IEEE frequency band.
        %
        %INPUT:
        % band: A string or character array which identifies the
        %       abbreviated IEEE frequency band.
        %
        %OUTPUT:
        % freq: A 1-by-2 vector containing the lower and upper bound for
        %       the specified frequency band specified in Hertz. If band
        %       does not match one of the known abbreviations, [NaN,NaN]
        %       is returned.
        %
        %December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        %
        bounds = [0.003,0.03,0.3,1,2,4,8,12,18,27,40,75,110,300]*1e9;
        bands = ["HF","VHF","UHF","L","S","C","X","Ku","K","Ka","V","W","mm"];
        i = find(bands==band);
        if any(i)
            freq = bounds(i:i+1);
        else
            freq = [NaN,NaN];
        end
    end
    function freq = getNATOBounds(band)
        %%GETNATOBOUNDS Returns the frequency bounds which define the given
        %               NATO frequency band.
        %
        %INPUT:
        % band: A string or character array which identifies the
        %       abbreviated NATO frequency band.
        %
        %OUTPUT:
        % freq: A 1-by-2 vector containing the lower and upper bound for
        %       the specified frequency band specified in Hertz. If band
        %       does not match one of the known abbreviations, [NaN,NaN]
        %       is returned.
        %
        %December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        %
        bounds = 1e6*[0,250,500,1e3,2e3,3e3,4e3,6e3,8e3,1e4,2e4,4e4,6e4,1e5,2e5,3e5];
        bands = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O"];
        i = find(bands==band);
        if any(i)
            freq = bounds(i:i+1);
        else
            freq = [NaN,NaN];
        end
    end
    function band = getITUBand(freq)
        %%GETITUBAND Returns the ITU frequency band in which the given
        %            frequency falls.
        %
        %INPUT:
        % freq: A scalar specifying the frequency in Hertz.
        %
        %OUTPUT:
        % band: A string which identifies the abbreviation of the
        %       corresponding ITU frequency band.
        %
        %December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        %
        bounds = [freq,0,3*10.^(1:13)];
        [~,I] = sort(bounds);
        bands = ["TLF","ELF","SLF","ULF","VLF","LF","MF","HF","VHF","UHF","SHF","EHF","THF"];
        if I(1)==1 || I(end)==1
            band = 'Undefined';
        else
            band = bands(I==1);
        end
    end
    function band = getIEEEBand(freq)
        %%GETIEEEBAND Returns the IEEE frequency band in which the given
        %             frequency falls.
        %
        %INPUT:
        % freq: A scalar specifying the frequency in Hertz.
        %
        %OUTPUT:
        % band: A string which identifies the abbreviation of the
        %       corresponding IEEE frequency band.
        %
        %December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        %
        bounds = [freq,[0.003,0.03,0.3,1,2,4,8,12,18,27,40,75,110,300]*1e9];
        [~,I] = sort(bounds);
        bands = ["HF","VHF","UHF","L","S","C","X","Ku","K","Ka","V","W","mm"];
        if I(1)==1 || I(end)==1
            band = 'Undefined';
        else
            band = bands(find(I==1)-1);
        end
    end
    function band = getNATOBand(freq)
        %%GETNATOBAND Returns the NATO frequency band in which the given
        %             frequency falls.
        %
        %INPUT:
        % freq: A scalar specifying the frequency in Hertz.
        %
        %OUTPUT:
        % band: A string which identifies the abbreviation of the
        %       corresponding NATO frequency band.
        %
        %December 2022 Codie T. Lewis, Naval Research Laboratory, Washington D.C.
        %
        bounds = [freq,1e6*[0,250,500,1e3,2e3,3e3,4e3,6e3,8e3,1e4,2e4,4e4,6e4,1e5]];
        [~,I] = sort(bounds);
        bands = ["A","B","C","D","E","F","G","H","I","J","K","L","M"];
        if I(1)==1 || I(end)==1
            band = 'Undefined';
        else
            band = bands(find(I==1)-1);
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
