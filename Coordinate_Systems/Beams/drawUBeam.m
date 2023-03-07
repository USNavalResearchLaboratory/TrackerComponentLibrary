function drawUBeam(uCenter,beamHalfwidth,range,lRx,M,closeEnds,varargin)
%%DRAWUBEAM Draw the outline in 2D Cartesian coordinates of a beam in a u
%           direction cosine in a particular direction. This does not shape
%           the beam based upon gain. It just outlines a region, drawing a
%           curve on the end at a constant range, if desired. This can be
%           useful for determining whether plotted targets are near various
%           beams in 2D. For beams that go outside of the -1, 1 unambiguous
%           region, their edges that would go outside the -1,1 region are
%           aliased back and two beams are drawn. The target is assumed to
%           be in front of the sensor, not behind it (positive y in the
%           local cooridnate system).
%
%INPUTS: uCenter The scalar center angle of the beam in direction cosine
%                space. This is from -1 to 1.
%  beamHalfwidth The distance from the center of the beam to one side of
%                the beam in direction cosine space. This must be positive
%                and values greater than or equal to 1 select the entire
%                visible region.
%          range This can be a scalar or a length-2 vector. If a scalar,
%                it is the maximum extent of the range to which the beam
%                is drawn. If this is a length-2 vector, then range(1) is
%                the minimum range to draw and range(2) is the maximum
%                range to be covered by the beam. This is a one-way range.
%            lRx The 2X1 [x;y] location vector of the receiver in
%                Cartesian coordinates.  If this parameter is omitted or
%                an empty matrix is passed, then the receiver is assumed
%                to be at the origin.
%              M A 2X2 rotation matrices to go from the alignment of the
%                global coordinate system to that at the receiver. If
%                omitted or an empty matrix is passed, then it is assumed
%                that the local coordinate system is aligned with the
%                global and M=eye(2) --the identity matrix is used.
%
%OUTPUTS: None.
%
%EXAMPLE:
%In this example, we have a beam on the boresight (red) and a beam of the
%same width that is so far from boresight that it ends up getting aliased
%back to the opposite direction. Note that in angular space, the beam of
%boresight is much wider than the beam that is centered on boresight. That
%is a property of using direction cosines.
% range=[2;10];
% lRx=[1;1];
% theta=deg2rad(30);
% M=rotMat2D(theta);
% beamHalfwidth=0.1;
% closeEnds=true;
% figure(1)
% clf
% hold on
% uCenter=0.0;
% drawUBeam(uCenter,beamHalfwidth,range,lRx,M,closeEnds,'-r','linewidth',2)
% uCenter=0.95;
% drawUBeam(uCenter,beamHalfwidth,range,lRx,M,closeEnds,'-k','linewidth',2)
% axis([-12,12,-12,12])
%
%January 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(lRx))
    lRx=[0;0];
end

if(nargin<5||isempty(M))
    M=eye(2,2);
end

if(nargin<6||isempty(closeEnds))
    closeEnds=true;
end

%Force it to be in the valid range.
uCenter=wrapRange(uCenter,-1,1);

if(beamHalfwidth<=0)
    error('beamHalfwidth must be positive.')
end

numPts=1000;

if(isscalar(range))
   %If only a maximum range is given, set the minimum range to 0.
   range=[0;range];
end

%Save the value of hold on the plot so that it can be reset to its previous
%state after plotting multiple lines.
holdVal=ishold();

if(beamHalfwidth>=1)
    %The special case where the entire viewing area is covered by the beam.
    uMin=-1;
    uMax=1;
else
    %First, check whether the region should be split into two regions. That
    %is, if the region goes below -1 in u, then it aliases over to the
    %positive side and if it goes above +1, then it aliases over to the
    %negative side.
    uMin=uCenter-beamHalfwidth;
    uMax=uCenter+beamHalfwidth;

    if(uMin<-1)
        %Split into two regions. One is from -1 to uMax. The other is from
        %the aliased uMin to 1.
        uMin=[-1,wrapRange(uMin,-1,1)];
        uMax=[uMax,1];
    elseif(uMax>1)
        %Split into two regions. One is from uMin to 1. The other is from
        %-1 to the aliased uMax.
        uMin=[uMin,-1];
        uMax=[1,wrapRange(uMax,-1,1)];
    end
end

numU=length(uMin);
useHalfRange=true;
%We are drawing either 1 or two regions, depending on whether or not the
%beam is aliased across the -1,1 boundary.
for curU=1:numU
    u1=uMin(curU);
    u2=uMax(curU);
    if(closeEnds==false)
        %Draw the sides
        startPoint=ru2Cart2D([range(1);u1],useHalfRange,lRx,lRx,M);
        endPoint=ru2Cart2D([range(2);u1],useHalfRange,lRx,lRx,M);
        plot([startPoint(1);endPoint(1)],[startPoint(2);endPoint(2)],varargin{:})
    
        hold on
        startPoint=ru2Cart2D([range(1);u2],useHalfRange,lRx,lRx,M);
        endPoint=ru2Cart2D([range(2);u2],useHalfRange,lRx,lRx,M);
        plot([startPoint(1);endPoint(1)],[startPoint(2);endPoint(2)],varargin{:})
        
        continue;
    end
    
    if(range(1)~=0)
        totalNumPoints=2*numPts+1;
    else
        totalNumPoints=numPts+2;
    end
    zPts=zeros(2,totalNumPoints);

    u=[u1+linspace(0,u2-u1,numPts-1),u2];
    zPts(:,1:numPts)=ru2Cart2D([range(2)*ones(1,numPts);u],useHalfRange,lRx,lRx,M);
    
    if(range(1)~=0)
        u=flip(u);
        zPts(:,(numPts+1):(2*numPts))=ru2Cart2D([range(1)*ones(1,numPts);u],useHalfRange,lRx,lRx,M);
    else
        zPts(:,numPts+1)=lRx;    
    end
    zPts(:,totalNumPoints)=zPts(:,1);   
    plot(zPts(1,:),zPts(2,:),varargin{:})
end

%Restore the hold value to its original setting.
if(~holdVal)
    hold off
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
