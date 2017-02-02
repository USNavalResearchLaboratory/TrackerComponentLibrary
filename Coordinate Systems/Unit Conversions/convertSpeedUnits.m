function val=convertSpeedUnits(val,unitOrigNum,unitOrigDenom,unitDestNum,unitDestDenom)
%%CONVERTSPEEDUNITS Convert values of speed from one set of units to
%                   another. Speed is given as a ratio of a distance unit
%                   to a time unit. Note that 1 knot (kt)=1 nml/h. Thus, a
%                   conversion from knots to meters per second would have 
%                   unitOrigNum='nml';unitOrigDenom='h';unitDestNum='m';
%                   unitDestDenom='s';
%
%INPUTS: val           The matrix or vector of values that are to be
%                      converted.
%unitOrigNum,unitOrigDenom
%unitDestNum,unitDestDenom Four character strings indicating the units of
%                      val and the units into which it is to be converted.
%                      unitOrigNum and unitDestNum are the units of
%                      distance and can take the values listed in
%                      convertTimeUnits. unitOrigDenom and unitDestDenom
%                      are the length units and can take the values listed
%                      in convertLengthUnits.
%
%OUTPUTS: val  The values converted into the desired coordinate system.
%
%This just uses the functions convertLengthUnits and convertTimeUnits to
%get the approprimate multiplication factors to convert the length and time
%components.
%
%May 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

lengthCoeff=convertLengthUnits(1,unitOrigNum,unitDestNum);
%Inverse, because time is in the denominator.
timeCoeff=1/convertTimeUnits(1,unitOrigDenom,unitDestDenom);

val=val*lengthCoeff*timeCoeff;
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
