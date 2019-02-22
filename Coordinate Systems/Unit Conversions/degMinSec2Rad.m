function rad=degMinSec2Rad(deg,min,sec)
%%DEGMINSEC2RAD Convert an angle given in degrees, arcminutes and
%               arcseconds to radians.
%
%INPUTS: deg A scalar or matrix number of degrees.
%        min A number of arcminutes, having the same dimensions as deg.
%        sec A number of arcseconds, having the same dimensions as deg.
%
%OUTPUTS: rad The angle or angles converted into radians, wrapped to -pi to
%             pi.
%
%There are 60 arcminutes per degree, 60 arcseconds per arcminute, and
%(180/pi) degrees per radians.
%
%October 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

deg=deg+min/60+sec/(60*60);
rad=wrapRange((pi/180)*deg,-pi,pi,false);
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
