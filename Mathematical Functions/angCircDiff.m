function diffVal=angCircDiff(ang)
%%ANGCIRCDIFF Given two angles (in radians) on a unit circle, take their
%         difference in radians. The difference can only span +/-pi. For
%         example the difference between 0 and 2*pi should be zero as both
%         represent the same point.
%
%INPUTS: angs A 2XnumAng matrix of angles in radians. The differences 
%               angs(1,:)-angs(2,:) in each column will be taken.
%
%OUTPUTS: diffVal A 1XnumAng vector of the differences on the unit circle
%                 of the angles. The shortest way around the circle is
%                 chosen. This goes from -pi to pi.
%
%This function simple takes the difference and wraps it to the range of -pi
%to pi with the wrapRange function.
%
%EXAMPLE:
%One can verify that for an angle of 0.1
% angCircDiff([0;0.1])-angCircDiff([0;0.1+2*pi])
%the difference is around machine precision. Similarly,
% angCircDiff([0;2*pi])
% is exactly zero
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

diffVal=wrapRange(ang(1,:)-ang(2,:),-pi,pi);

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
