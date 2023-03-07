function plotBistaticRangeEllipse(rMeas,sigmaR2,lRx,lTx,lRef,MRef,refPoint,refType,gammaVal,patchOpts,varargin)
%%PLOTBISTATICRANGEELLIPSE Given a bistatic range measurement, plot a slice
%           of the bistatic range ellipse with respect to a particular
%           reference height either given explicitly or derived from a
%           reference point in a particularly rotated coordinate system.
%           This can be useful for visualizaing trilateration problems. One
%           can also pass a standard deviation of the range measurement in
%           which case two ellipses defining a confidence interval are
%           plotted.
%
%INPUTS: rMeas The scalar bistatic range measurement.
%      sigmaR2 The variance of the bistatic range measurement. If an empty
%              matrix is passed here, then only a line at the measured
%              bistatic range is drawn rather than two lines marking a
%              confidence region.
%          lRx The 3X1 location of the receiver in global coordinates.
%          lTx The 3X1 location of the transmitter in global coordinates.
%    lRef,MRef A 3X1 point and a 3X3 rotation matrix that define how the
%              coordinate system used for plotting the measurements differs
%              from the global cooridnate system. lRef defined the local
%              origin and it is the point about which a rotation by MRef
%              will rotate from the global coordinate system into the local
%              coordinate system for plotting. This rotation can be
%              important, because the plane used to cut the 3D ellipsoid to
%              make a 2D plot is the x-y plane (at a particular z height)
%              in the _local_ coordinate system. If omitted or empty
%              matrices are passed, lRef=lRx and MRef=eye(3,3).
%     refPoint This can either be a scalar value or a 3X1 point. If a
%              scalar, value, this defines the height in the local z axis
%              at which the x-y plane is cut through the bistatic range
%              ellipsoid. If this is a 3X1 vector, then if refType=0, the
%              height (z value) of the point transformed into the local
%              coordinate system is used to define where the x-y cutting
%              plane is. On the other hand, if refType=1, then for each
%              ellipsoid, the point on the ellipsoid that is closest to
%              refPoint is used as the reference point. If omitted or an
%              empty matrix is passed, then refPoint=lRx is used.
%      refType This specifies how refPoint is used, as described above. If
%              omitted or an empty matrix is passed, then refType=0 is
%              used. refType must be 0 if only a scalar height is passed
%              for refPoint.
%     gammaVal This is only used if sigmaR2 is passed. This is the
%              threshold defining the confidence region to plot. The
%              confidence region is defined as
%              (r-rMeas)^2/sigmaR2<=gammaVal. Thus, the minimum and maximum
%              ranges that will be plotted are
%              r=rMeas+/-sqrt(gammaVal*sigmaR2)
%              If omitted or an empty matrix is passed, the value
%              corresponding to the 99.97% confidence region is used and
%              obtained as ChiSquareD.invCDF(0.9997,1).
%    patchOpts If sigmaR2 is provided,then this input is used (otherwise an
%              empty matrix can be passed or it is omitted). If this input
%              is -1, then the space between the confidence ellipses is not
%              colored in. If one wishes to color in the values between the
%              ellipses, then if this parameter is provided, then this is a
%              cell array. patchOpts{1} must be an RGB color that the space
%              in between will be colored in with. The other elements are
%              then options for the Matlab patch function. The default if
%              omitted or an empty matrix is passed is
%              {[0.5,0.5,0.5],'linestyle','non','facealph',0.5};
%     varargin This means that one can add on any number of comma-separated
%              options that will be passed as inputs to the drawEllips
%              function for plotting.
%
%OUTPUTS: None.
%
%The function puts the bistatic ellipse parameters into centered quadratic
%form using getBistaticEllipse. The cut with the x-y plane is performed
%using the projEllipse2ZPlane function. Plotting is then done using the 
%drawEllipse function.
%
%EXAMPLE 1:
%This is an example of showing how multiple bistatic ellipses intersect at
%the location of a target in the noise-free case.
% %The transmitter locations.
% xTx=[[0;5e3;8e3],[0;0;10e3],[0;0;0],[0;0;5e4]];
% %The receiver locations.
% xRx=[[0;0;0],[0;2e4;0],[2e4;0;0],[0;-2e4;0]];
% S=size(xRx,2);%The number of bistatic paths.
% %The target location.
% t=[-18500;10000;8000];
% %Generate range measurements. In this example, we do not add noise.
% rMeas=zeros(1,S);
% useHalfRange=false;
% for k=1:S
%     rMeas(1,k)=getRange(t,useHalfRange,xTx(:,k),xRx(:,k));
% end
% 
% %Plot the bistatic ellipses for the given range and transmitter and
% %receiver locations. We plot then all projected onto the plane at the
% %altitude of the target.
% figure(1)
% clf
% hold on
% colors={'-c','-m','-r','-k'};
% for k=1:S
%     plotBistaticRangeEllipse(rMeas(1,k),[],xRx(:,k),xTx(:,k),xRx(:,1),eye(3,3),t,[],[],colors{k},'linewidth',2)
% end
% scatter(t(1),t(2),400,'.g')%Mark the target location.
%
%EXAMPLE 2:
%This is similiar to the first example, except noise is added to the
%measurements and their confidence regions are plotted rather than just
%drawing lines at the measurement locations. In this instance, we color in
%the space between the lines.
% %The transmitter locations.
% xTx=[[0;5e3;8e3],[0;0;10e3],[0;0;0],[0;0;5e4]];
% %The receiver locations.
% xRx=[[0;0;0],[0;2e4;0],[2e4;0;0],[0;-2e4;0]];
% S=size(xRx,2);%The number of bistatic paths.
% %The target location.
% t=[-18500;10000;8000];
% sigmaR=1e3;
% %Generate range measurements.
% rMeas=zeros(1,S);
% useHalfRange=false;
% for k=1:S
%     rMeas(1,k)=getRange(t,useHalfRange,xTx(:,k),xRx(:,k))+sigmaR*randn(1);
% end
% %Plot the bistatic ellipses for the given range and transmitter and
% %receiver locations. We plot then all projected onto the plane at the
% %altitude of the target.
% figure(1)
% clf
% hold on
% colors={'c','m','r','k'};
% patchOpts{1}={[0,1,1],'linestyle','non','facealph',0.5};
% patchOpts{2}={[1,0,1],'linestyle','non','facealph',0.5};
% patchOpts{3}={[1,0,0],'linestyle','non','facealph',0.5};
% patchOpts{4}={[0,0,0],'linestyle','non','facealph',0.5};
% 
% for k=1:S
%     plotBistaticRangeEllipse(rMeas(1,k),sigmaR^2,xRx(:,k),xTx(:,k),xRx(:,1),eye(3,3),t,[],[],patchOpts{k},colors{k},'linewidth',2)
% end
% scatter(t(1),t(2),400,'.g')%Mark the target location.
%
%December 2020 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<10||isempty(patchOpts))
   patchOpts={[0.5,0.5,0.5],'linestyle','non','facealph',0.5}; 
end

if(nargin<9||isempty(gammaVal))
    gammaVal=ChiSquareD.invCDF(0.9997,1);
end

if(nargin<8||isempty(refType))
    refType=0; 
end

if(nargin<7||isempty(refPoint))
    refPoint=lRx;
end

if(nargin<6||isempty(MRef))
    MRef=eye(3,3);
end

if(nargin<5||isempty(lRef))
    lRef=lRx;
end

numPoints=2000;

%Get the maximum and minimum bistatic ranges.
if(isempty(sigmaR2))
    r=rMeas;
    numR=1;
else
    deltaR=sqrt(gammaVal*sigmaR2);
    r=zeros(2,1);
    r(1)=rMeas+deltaR;%Maximum radius.
    r(2)=rMeas-deltaR;%Minimum radius.
    numR=2;
    zPoints=zeros(2,numPoints,2);
end

%Put everything into the receiver's local coordinate system.
lRx=MRef*(lRx-lRef);
lTx=MRef*(lTx-lRef);
if(isscalar(refPoint))
    h=refPoint;
else
    refPoint=MRef*(refPoint-lRef);
    h=refPoint(3);%The height to use for the 2D slice.
end

%Save the value of hold on the plot so that it can be reset to its previous
%state after plotting possibly two ellipses.
holdVal=ishold();
hold on
for k=1:numR
    %Put the bistatic ellipsoid into centered quadratic form.
    [M3D,kpp3D,cNew3D]=getBistaticEllipse(r(k),lTx,lRx);
    
    if(refType==1)
        %If the height should be set by the point on the outer ellipsoid
        %that is closest to the point of interest.
        zp=nearestPointOnEllipsoid(kpp3D,M3D,refPoint,cNew3D);
        h=zp(3);
    end

    %Project onto an x-p plane at the true height of the target.
    [M,kpp,c]=projEllipse2ZPlane(M3D,kpp3D,cNew3D,h);
    
    if(~isempty(M))
        %If the ellipsoid intersects the plane.
        if(~isempty(sigmaR2)&&patchOpts~=-1)
            zPoints(:,:,k)=getEllipsePoints(kpp,M,c,numPoints,false);
        end

        drawEllipse(kpp,M,c,varargin{:})
    else
        %Restore the hold value to its original setting.
        if(~holdVal)
            hold off
        end
        return;
    end
end

if(~isempty(sigmaR2)&&patchOpts~=-1)
    x=zPoints(1,:,1);
    X=zPoints(1,:,2);
    y=zPoints(2,:,1);
    Y=zPoints(2,:,2);

    patch([x X],[y Y],patchOpts{:});
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
