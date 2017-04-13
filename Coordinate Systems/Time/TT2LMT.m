function [LMTRad,LMT1,LMT2]=TT2LMT(TT1,TT2,rObsITRS,deltaT,xpyp)
%%TT2LMT Convert from terrestrial time to local mean solar time (LMT). This
%        is a measure of the time at a place on the Earth given the mean
%        location of the Sun, not a high-precision ephemeris.
%
%INPUTS: Jul1,Jul2 Two parts of a Julian date given in TT. The
%                  units of the date are days. The full date is the sum of
%                  both terms. The date is broken into two parts to
%                  provide more bits of precision. It does not matter how
%                  the date is split.
%         rObsITRS The 3X1 location of the observer in the International
%                  Terrestrial Reference System (ITRS). Only the direction
%                  matters, not the magnitude.
%           deltaT An optional parameter specifying the offset between TT
%                  and UT1 in seconds. If this parameter is omitted or if
%                  an empty matrix is passed, then the value of the
%                  function getEOP will be used.
%             xpyp xpyp=[xp;yp] are the polar motion coordinates in radians
%                  including the effects of tides and librations. If this
%                  parameter is omitted or if an empty matrix is passed,
%                  the value from the function getEOP will be used.
%
%OUTPUTS: LMTRad The local apparent solar time in radians. 0 radians is
%                solar midnight. pi radians is solar noon. The mapping is
%                is 2*pi radians for 24 hours. Thus, LATRad*(24/(2*pi))
%                gives the time of day in hours.
%      LMT1,LMT2 Two parts of the local mean solar time in Julian days. The
%                date is split so that LMT1 is the integer part and LMT2 is
%                the fractional part. Like UT1, a zero fractional part of
%                the day corresponds to noon, not midnight.
%
%The local mean solar time is UT1 plus the East longitude of the observer
%in radians in the Terrestrial Intermediate Reference System (TIRS). The
%East longitude is added using the convention that 360 degrees equals 24
%hours of 60 minutes and 60 seconds. For LMT in radians, an additional 12
%hours (pi radians) is added to adhere to the convention that 0 hours/ 0
%radians is midnight, not noon.
%
%The use of UT1 as a definition of Greenwhich mean solar time is given in
%[1]. Since UT1 is defined as a rotation in the TIRS, it only makes sense
%that a local mean solar time is defined by using the longitude offset in
%the TIRS.
%
%REFERENCES:
%[1] D. D. McCarthy, "Evolution of timescales from astronomy to physical
%    metrology," Metrologia, vol. 48, no. 4, pp. S132-S144, Aug. 2011.
%
%April 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Get any Earth orientation parameters that were not provided.
if(nargin<5||isempty(deltaT)||isempty(xpyp))
    [UTC1,UTC2]=TT2UTC(TT1,TT2);
    [xpypNew,~,~,deltaTTUT1]=getEOP(UTC1,UTC2);
    
    if(nargin<4||isempty(deltaT))
       deltaT=deltaTTUT1;
    end
    
    if(nargin<5||isempty(xpyp))
        xpyp=xpypNew;
    end
end

[UT11,UT12]=TT2UT1(TT1,TT2,deltaT);
rObsTIRS=ITRS2TIRS(rObsITRS,TT1,TT2,xpyp);
rObsSphere=Cart2Sphere(rObsTIRS);

%2*pi radians per day.
deltaUT1=rObsSphere(2)/(2*pi);

%Add preserving precision.
if(UT11<UT12)
    LMT2=UT11+deltaUT1;
    LMT1=UT12;
else
    LMT2=UT12+deltaUT1;
    LMT1=UT11;
end

%Put the fractional part into LMT2 and the integer part into LMT1.
frac1=LMT1-fix(LMT1);
frac2=LMT2-fix(LMT2);
sumVal=frac1+frac2;
sumInt=fix(sumVal);
sumFrac=sumVal-sumInt;
LMT1=fix(LMT2)+fix(LMT1)+sumInt;
LMT2=sumFrac;

%The local mean solar time in radians. This is equivalent to
%wrapRange(2*pi*(LMT1+LMT2)+pi,-pi,pi); However, since LMT1 is always
%integers, it can just be omitted.
LMTRad=wrapRange(2*pi*LMT2+pi,-pi,pi);

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
