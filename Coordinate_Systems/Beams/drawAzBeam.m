function drawAzBeam(azCenter,beamHalfwidth,range,systemType,lRx,M,closeEnds,varargin)
%%DRAWANGLEBEAM Draw the outline in 2D Cartesian coordinates of a polar
%           beam in a particular direction. This does not shape the beam
%           based upon gain. It just outlines a region, drawing a curve on
%           the end at a constant range, if desired. This can be useful for
%           determining whether plotted targets are near various beams in
%           2D.
%
%INPUTS: azCenter The scalar center angle of the beam in radians.
%   beamHalfwidth The distance from the center of the beam to one side of
%                 the beam in radians. This must be >0.
%           range This can be a scalar or a length-2 vector. If a scalar,
%                 it is the maximum extent of the range to which the beam
%                 is drawn. If this is a length-2 vector, then range(1) is
%                 the minimum range to draw and range(2) is the maximum
%                 range to be covered by the beam. This is a one-way range.
%      systemType An optional parameter specifying the axis from which the
%                 angles are measured. Possible values are:
%                 0 (The default if omitted or an empty matrix is passed)
%                   The azimuth angle is counterclockwise from the x axis.
%                 1 The azimuth angle is measured clockwise from the y axis.
%             lRx The 2X1 [x;y] location vector of the receiver in
%                 Cartesian coordinates.  If this parameter is omitted or
%                 an empty matrix is passed, then the receiver is assumed
%                 to be at the origin.
%               M A 2X2 rotation matrices to go from the alignment of the
%                 global coordinate system to that at the receiver. If
%                 omitted or an empty matrix is passed, then it is assumed
%                 that the local coordinate system is aligned with the
%                 global and M=eye(2) --the identity matrix is used.
%
%OUTPUTS: None.
%
%EXAMPLE:
%This function draws 60-degree beams in three directions making an image
%reminiscent of a radioactive symbol.
% figure(1)
% clf
% hold on
% %60 degree beam.
% beamHalfwidth=30*(pi/180);
% range=[2;10];
% systemType=0;
% azCenter=30*(pi/180);
% drawAzBeam(azCenter,beamHalfwidth,range,[],[],[],[],'-k','linewidth',2)
% azCenter=pi-30*(pi/180);
% drawAzBeam(azCenter,beamHalfwidth,range,[],[],[],[],'-k','linewidth',2)
% azCenter=-pi/2;
% drawAzBeam(azCenter,beamHalfwidth,range,[],[],[],[],'-k','linewidth',2)
%
%January 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<4||isempty(systemType))
    systemType=0; 
end

if(nargin<5||isempty(lRx))
    lRx=[0;0];
end

if(nargin<6||isempty(M))
    M=eye(2,2);
end

if(nargin<7||isempty(closeEnds))
    closeEnds=true;
end

numPts=1000;

if(isscalar(range))
   %If only a maximum range is given, set the minimum range to 0.
   range=[0;range];
end

if(systemType==1)
    azCenter=pi/2-azCenter;
end

if(beamHalfwidth<=0)
    error('beamHalfwidth must be positive.')
end

%Save the value of hold on the plot so that it can be reset to its previous
%state after plotting multiple lines.
holdVal=ishold();

systemType=0;
useHalfRange=true;

if(closeEnds==false)
    %Draw the sides
    az1=azCenter-beamHalfwidth;
    startPoint=pol2Cart([range(1);az1],systemType,useHalfRange,lRx,lRx,M);
    endPoint=pol2Cart([range(2);az1],systemType,useHalfRange,lRx,lRx,M);
    plot([startPoint(1);endPoint(1)],[startPoint(2);endPoint(2)],varargin{:})

    hold on
    az2=azCenter+beamHalfwidth;
    startPoint=pol2Cart([range(1);az2],systemType,useHalfRange,lRx,lRx,M);
    endPoint=pol2Cart([range(2);az2],systemType,useHalfRange,lRx,lRx,M);
    plot([startPoint(1);endPoint(1)],[startPoint(2);endPoint(2)],varargin{:})
    
    %Restore the hold value.
    if(~holdVal)
        hold off
    end
    return
end

if(range(1)~=0)
    totalNumPoints=2*numPts+1;
else
    totalNumPoints=numPts+2;
end
zPts=zeros(2,totalNumPoints);
    
az1=azCenter-beamHalfwidth;
az2=azCenter+beamHalfwidth;
ang=[az1+linspace(0,2*beamHalfwidth,numPts-1),az2];
zPts(:,1:numPts)=pol2Cart([range(2)*ones(1,numPts);ang],systemType,useHalfRange,lRx,lRx,M);

if(range(1)~=0)
    ang=flip(ang);
    zPts(:,(numPts+1):(2*numPts))=pol2Cart([range(1)*ones(1,numPts);ang],systemType,useHalfRange,lRx,lRx,M);
else
    zPts(:,numPts+1)=lRx;    
end
zPts(:,totalNumPoints)=zPts(:,1);   
plot(zPts(1,:),zPts(2,:),varargin{:})

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
