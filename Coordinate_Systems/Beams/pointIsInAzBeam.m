function boolVal=pointIsInAzBeam(points,azBounds,boundType)
%%POINTISINAZBEAM A region in angles is defined either in
%      terms of a central angle and a beamwidth (the width being the
%      distance between the angular bounds of the region) or in terms of a
%      starting angle and an ending angle, determine whether provided
%      points are within that angular span. Because of the circular nature
%      of angles, one cannot simply compare the values in the points to
%      those of the bounds. 
%
%INPUTS: points A vector or matrix of angular values in radians to test
%               whether or not they are in the beam.
%      azBounds A 2X1 vector defining the bounds of the beam. If
%               boundType=0, then azBounds(1) is the direction of the
%               center of the beam in radians and azBounds(2) is the half
%               beamwidth of the beam. If boundType=1, then azBounds(1) is
%               the lower angular bound and azBounds(2) is the upper
%               angular bound --where the angles are taken in increasing
%               (aliased) order. Thus, if azBounds(1) is 0 and azBounds(2)
%               is -pi/2, then the beam occupies 3/4 of the circle.
%     boundType This is either 0 or 1 and affects the format of azBounds,
%               as described above. The default if omitted or an empty
%               matrix is passed is 0. 
%
%OUTPUTS: boolVal This has the same size as the input and the entries are
%                 true if the corresponding entry in the input is in the
%                 beam and false if it is not in the beam. 
%
%EXAMPLE:
%This takes all the angles on the unit circle and plots them. It also forms
%a beam and plots just the points in the beam. Points are determined to be
%in the beam or not using pointIsInAzBeam.
% beamCenter=deg2rad(130);
% beamHalfwidth=deg2rad(15);
% numPts=1000;
% points=linspace(-pi,pi,numPts);
% boolVals=pointIsInAzBeam(points,[beamCenter;beamHalfwidth]);
% 
% systemType=0;
% xyPts=pol2Cart(points,systemType);
% figure(1)
% clf
% hold on
% closeEnds=true;
% drawAzBeam(beamCenter,beamHalfwidth,1,systemType,[],[],closeEnds,'-k','linewidth',4);
% scatter(xyPts(1,:),xyPts(2,:),'.b')
% axis([-1.25,1.25,-1.25,1.25])
% %Only the points in the beam.
% xyPts=xyPts(:,boolVals);
% scatter(xyPts(1,:),xyPts(2,:),'.r')
%
%January 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(boundType))
    boundType=0;
end

if(boundType==0)
    %If defined in terms of a beam center and a beamwidth....
    beamCenter=azBounds(1);
    beamHalfwidth=azBounds(2);

    if(beamHalfwidth>=pi)
        boolVal=true(size(points));
        return;
    end

    az1=beamCenter-beamHalfwidth;
    az2=beamCenter+beamHalfwidth;
else
    %If the bounds are directly given.
    az1=azBounds(1);
    az2=azBounds(2);
end
%Rotate everything so that az1=0. Everything is considered in terms of
%increasing angle from 0.
az2=wrapRange(az2-az1,0,2*pi);
points=wrapRange(points-az1,0,2*pi);

%Test whether the angles are in or on the edge of the beam.
boolVal=(points<=az2);

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
