function boolVal=pointIsInUBeam(points,uBounds,boundType)
%%POINTISINUBEAM A region in a single direction cosine in 2D space is
%      defined either in terms of a central direction and a half beamwidth
%      (the beamwidth being the distance between the bounds of the region)
%      or in terms of a starting direction cosine and an ending direction
%      cosine, determine whether provided points are within that span.
%      Because beams and points that go past the -1 to 1 region alias back
%      into that region, one cannot simply compare the values in the points
%      to those of the bounds. 
%
%INPUTS: points A vector or matrix of direction cosines values to test
%               whether or not they are in the beam. nonaliased values are
%               between -1 and 1.
%       uBounds A 2X1 vector defining the bounds of the beam. If
%               boundType=0, then uBounds(1) is the direction of the
%               center of the beam and uBounds(2) is the half
%               beamwidth of the beam. If boundType=1, then uBounds(1) is
%               the lower bound and azBounds(2) is the upper bound --where
%               the values are taken in increasing (aliased) order. If
%               bounds go outside of the beam span, then they are assumed
%               to alias back into the valid span (a beam split between two
%               directions 180 degrees apart).
%     boundType This is either 0 or 1 and affects the format of uBounds,
%               as described above. The default if omitted or an empty
%               matrix is passed is 0. 
%
%OUTPUTS: boolVal This has the same size as the input and the entries are
%                 true if the corresponding entry in the input is in the
%                 beam and false if it is not in the beam. 
%
%EXAMPLE:
%This plots a beam that is aliased across the -1,1 boundary. Then, points
%uniformly spaced across the -1 to 1 u span are tested for being in the
%beam. Those in the beam are plotted in red and all the other points are
%plotted in blue. One can see that the function correctly handles aliasing
%around the -1,1 discontinuity.
% range=[0;1];
% beamHalfwidth=0.1;
% closeEnds=true;
% figure(1)
% clf
% hold on
% uCenter=0.95;
% drawUBeam(uCenter,beamHalfwidth,range,[],[],closeEnds,'-k','linewidth',2)
% 
% numPoints=1000;
% uPoints=linspace(-1,1,numPoints);
% boundType=0;
% uBounds=[uCenter;beamHalfwidth];
% boolVal=pointIsInUBeam(uPoints,uBounds,boundType);
% 
% %Plot all the points at unit range in blue.
% points=[ones(1,numPoints);uPoints];
% useHalfRange=true;
% xyPts=ru2Cart2D(points,useHalfRange);
% scatter(xyPts(1,:),xyPts(2,:),100,'.b')
% %Plot the points that gate in red.
% points=points(:,boolVal);
% xyPts=ru2Cart2D(points,useHalfRange);
% scatter(xyPts(1,:),xyPts(2,:),200,'.r')
% axis([-1,1,-0.5,1])
%
%January 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(boundType))
    boundType=0;
end

if(boundType==0)
    %If defined in terms of a beam center and a beamwidth....
    beamCenter=uBounds(1);
    beamHalfwidth=uBounds(2);

    if(beamHalfwidth>=1)
        boolVal=true(size(points));
        return;
    end

    u1=beamCenter-beamHalfwidth;
    u2=beamCenter+beamHalfwidth;
else
    %If the bounds are directly given.
    u1=uBounds(1);
    u2=uBounds(2);
end
%Alias everything so that u1=0. Everything is considered in terms of
%increasing u from 0.
u2=wrapRange(u2-u1,0,2);
points=wrapRange(points-u1,0,2);

%Test whether the direction cosines are in or on the edge of the beam.
boolVal=(points<=u2);

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
