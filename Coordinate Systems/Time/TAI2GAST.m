function GMST=TAI2GAST(Jul1,Jul2,version,deltaT)
%%TAI2GAST Convert from international atomic time (TAI) to Greenwhich
%          apparent sidereal time (GAST), which is a measure of the
%          rotational direction of the Earth.
%
%INPUTS: Jul1,Jul2 Two parts of a pseudo-Julian date given in TAI. The
%                  units of the date are days. The full date is the sum of
%                  both terms. The date is broken into two parts to
%                  provide more bits of precision. It does not matter how
%                  the date is split.
%          version An optional integer specifying the theory to use for
%                  GMST. The theory chosen should be consistent with other
%                  values used in astronomical routines. Possible values
%                  are
%                  1982 Compute GAST ion accordance with the International
%                     Astronomical Union's (IAU's) 1982 model.
%                  2000 Compute GMST in line with IAU 2000 resolutions
%                     related to precession and nutation.
%                  2006 (The default if omitted) Compute GMST in line with
%                     IAU 2006 resolutions related to precession and
%                     nutation.
%           deltaT An optional parameter specifying the offset between TT
%                  and UT1 in seconds. If this parameter is omitted, then
%                  the value of the function getEOP will be used.
%
%OUTPUTS: GAST The Greenwhich apparent sideral time in radians. 
%
%The function just calls TAI2TT and then TT2GAST.
%
%GAST is defined in Section 5.5.7 of [1].
%
%REFERENCES:
%[1] G. Petit and B. Luzum, IERS Conventions (2010), International Earth
%    Rotation and Reference Systems Service Std. 36, 2010.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

[Jul1,Jul2]=TAI2TT(Jul1,Jul2);

if(nargin==2)
    GMST=TT2GAST(Jul1,Jul2);
elseif(nargin==3)
    GMST=TT2GAST(Jul1,Jul2,version);
else
    GMST=TT2GAST(Jul1,Jul2,version,deltaT);
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
