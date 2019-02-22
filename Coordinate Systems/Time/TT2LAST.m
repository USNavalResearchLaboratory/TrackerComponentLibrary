function LAST=TT2LAST(Jul1,Jul2,rObsITRS,version,deltaT,xpyp)
%%TT2LAST Convert from terrestrial time (TT) to local apparent sidereal
%         time (LAST), which is a measure of the rotational direction of
%         the Earth.
%
%INPUTS: Jul1,Jul2 Two parts of a Julian date given in TT. The
%                  units of the date are days. The full date is the sum of
%                  both terms. The date is broken into two parts to
%                  provide more bits of precision. It does not matter how
%                  the date is split.
%         rObsITRS The 3X1 location of the observer in the International
%                  Terrestrial Reference System (ITRS). Only the direction
%                  matters, not the magnitude.
%          version An optional integer specifying the theory to use for
%                  GMST. The theory chosen should be consistent with other
%                  values used in astronomical routines. Possible values
%                  are
%                  1982 Compute GMST ion accordance with the International
%                     Astronomical Union's (IAU's) 1982 model.
%                  2000 Compute GMST in line with IAU 2000 resolutions
%                     related to precession and nutation.
%                  2006 (The default if omitted) Compute GMST in line with
%                     IAU 2006 resolutions related to precession and
%                     nutation.
%           deltaT An optional parameter specifying the offset between TT
%                  and UT1 in seconds. If this parameter is omitted or if
%                  an empty matrix is passed, then the value of the
%                  function getEOP will be used.
%             xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
%                  including the effects of tides and librations. If this
%                  parameter is omitted or if an empty matrix is passed,
%                  the value from the function getEOP will be used.
%
%OUTPUTS: LAST The local apparent sidereal time in radians.
%
%The relation between LMST and Greenwhich apparent sidereal time (GAST) is
%documented in [1]. LAST is GAST added to the East longitude of the
%observer in the TIRS.
%
%REFERENCES:
%[1] G. H. Kaplan, "The IAU resolutions on astronomical reference systems,
%    time scales, and Earth rotation models: Explanation and
%    implementation," U.S. Naval Observatory, Tech. Rep. 179, 20 Oct. 2005.
%    [Online].
%    Available: http://aa:usno:navy:mil/publications/docs/Circular 179:pdf
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4)
    version=2006;
end

%Get any Earth orientation parameters that were not provided.
if(nargin<6||isempty(deltaT)||isempty(xpyp))
    [UTC1,UTC2]=TT2UTC(Jul1,Jul2);
    [xpypNew,~,~,deltaTTUT1]=getEOP(UTC1,UTC2);
    
    if(nargin<5||isempty(deltaT))
       deltaT=deltaTTUT1;
    end
    
    if(nargin<6||isempty(xpyp))
        xpyp=xpypNew;
    end
end

%Get GAST.
GAST=TT2GAST(Jul1,Jul2,version,deltaT);
rObsTIRS=ITRS2TIRS(rObsITRS,Jul1,Jul2,xpyp);
rObsSphere=Cart2Sphere(rObsTIRS);

%Add the East longitude of the observer in the TIRS.
LAST=wrapRange(GAST+rObsSphere(2),-pi,pi);

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
