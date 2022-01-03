function boolVals=pointIsInPolAngSpan(angVals,angMin,angMax)
%%POINTISINPOLANGSPAN Given that an angular region goes from angMin to
%               angMax, determine if a particular direction is in the
%               region. Numerically, angMin can be less than angMax, so
%               there is no issue with cross a -pi/pi boundary. However,
%               the span is taken to be in the direction of increasing
%               angles from angMin. Thus, it is possible for the region to
%               be more than pi radians wide.
%
%INPUTS: angVals A matrix of angles to determine whether they are in the
%               beam.
% angMin, angMax The starting and ending angles of the beam. These must be
%               in the range of -pi to pi radians.
%
%OUTPUTS: boolVals A matrix having the same dimensions as angVals that
%                  indicates whether (true) or not (false) each value in
%                  angVals is in the beam.
%
%July 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

diff=angMax-angMin;
angSpan=wrapRange(diff,0,2*pi);
if(angSpan==0&&diff~=0)
    %It must go around at least once.
    angSpan=2*pi;
end

angVals=wrapRange(angVals-angMin,0,2*pi);
boolVals=angVals<=angSpan;

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
