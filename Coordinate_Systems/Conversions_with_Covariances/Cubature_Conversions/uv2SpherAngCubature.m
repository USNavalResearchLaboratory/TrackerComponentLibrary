function [azEl,RAzEl]=uv2SpherAngCubature(zUV,SR,systemType,Ms,Muv,xi,w)
%%UV2SPHERANGCUBATURE Approximate the mean and covariance matrix of a
%     Gaussian noise-corrupted direction measurement in u-v direction
%     cosine coordinates converted into spherical coordinates using
%     cubature integration.
%
%INPUTS: zUV A 2XnumMeas (for just u-v coordinate) set of N direction
%            cosines to convert.
%          R The 2X2XnumMeas measurement covariance matrices for the
%            measurements. If all of the matrices are the same, then this
%            can just be a single 2X2 matrix.
% systemType An optional parameter specifying the axis from which the
%           angles are measured in radians. Possible values are
%           0 (The default if omitted) Azimuth is measured 
%             counterclockwise from the x-axis in the x-y plane. Elevation
%             is measured up from the x-y plane (towards the z-axis). This
%             is consistent with common spherical coordinate systems for
%             specifying longitude (azimuth) and geocentric latitude
%             (elevation).
%           1 Azimuth is measured counterclockwise from the z-axis in the
%             z-x plane. Elevation is measured up from the z-x plane
%             (towards the y-axis). This is consistent with some spherical
%             coordinate systems that use the z-axis as the boresight
%             direction of the radar.
%           2 This is the same as 0 except instead of being given
%             elevation, one is given the angle away from the z-axis, which
%             is (pi/2-elevation).
%           3 This is the same as 0 except azimuth is measured clockwise
%             from the y-axis in the x-y plane instead of counterclockwise
%             from the x-axis. This coordinate system often arises when
%             given "bearings" in a local East-North-Up coordinate system,
%             where the bearing directions are measured East of North.
%   Ms, Muv If either the spherical coordinate system or the u-v coordinate
%           system is rotated compared to the global Cartesian coordinate
%           system, these optional 3X3 matrices provide the rotations. Ms
%           is a 3X3 matrix to go from the alignment of a global
%           Cartesian coordinate system to that in which the spherical
%           coordinates are computed. Similarly, Muv is a rotation matrix
%           to go from the alignment of a global Cartesian cordinate system
%           to that in which the u-v(-w) coordinates are computed. If
%           either of these in omitted or an empty matrix is passed, then
%           the missing one is replaced with the identity matrix.
%        xi A 2XnumCubaturePoints matrix of cubature points for the numeric
%           integration. If this and the final parameter are omitted or
%           empty matrices are passed, then fifthOrderCubPoints is used to
%           generate cubature points.
%         w A numCubaturePointsX1 vector of the weights associated with the
%           cubature points.
% 
%OUTPUTS: zCart The approximate means of the PDF of the spherical converted
%               measurements in [azimuth;elevation] coordinates in radians
%               for each measurement. This is a 2XnumMeas matrix.
%         RCart The approximate 2X2XnumMeas set of covariance matrices of
%               the PDFs of the spherical converted measurements. 
%
%The general idea behind cubature integration measurement conversion is
%described in [1]. The azimuth is wrapped to avoid issues around +/- pi
%boundaries. However, the elevation is not wrapped, so it is assumed that
%the target is not near the extremeties in elevation.
%
%EXAMPLE:
%In this examples the NEES of converted measurents with a fixed v and a
%changing u is shown. The NEES is shown to be within 99.9% confidence
%bounds.
% numMCRuns=1e3;
% numUVals=50;
% systemType=0;
% v=0.4;%A fixed v value.
% uMax=sqrt(1-v^2);
% uVals=linspace(-0.6*uMax,0.6*uMax,numUVals);
% R=diag([2e-3;2e-3]);
% SR=chol(R,'lower');
% Ms=randRotMat(3);
% Muv=Ms;
% uvVals=[uVals;v*ones(1,numUVals)];
% [xi,w]=fifthOrderCubPoints(2);
% NEES=zeros(numUVals,1);
% for k=1:numUVals
%     zCurTrue=uvVals(:,k);
%     zAzElTrue=uv2SpherAng(zCurTrue,systemType,Ms,Muv);
% 
%     for curRun=1:numMCRuns
%         z=zCurTrue+SR*randn(2,1);
%         [zAzEl,RAzEl]=uv2SpherAngCubature(z,SR,systemType,Ms,Muv,xi,w);
%         diffVal=zAzEl-zAzElTrue;%Ignoring wrapping.
%         NEES(k)=NEES(k)+diffVal'*inv(RAzEl)*diffVal;
%     end
% end
% NEES=NEES/(2*numMCRuns);
% bounds=getNEESConfBounds(0.999,2,numMCRuns);
% figure(1)
% clf
% hold on
% plot(uVals,NEES,'linewidth',2)
% %Plot the lines for the 99.9% confidence region.
% plot([uVals(1),uVals(end)],[bounds(1),bounds(1)],'-k','linewidth',2)
% plot([uVals(1),uVals(end)],[bounds(2),bounds(2)],'-k','linewidth',2)
% axis([uVals(1),uVals(end),0,2])
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%May 2024 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<6||isempty(xi))
    [xi,w]=fifthOrderCubPoints(2);
end

if(nargin<5||isempty(Muv))
    Muv=eye(3,3);
end

if(nargin<4||isempty(Ms))
    Ms=eye(3,3);
end

if(nargin<3||isempty(systemType))
    systemType=0;
end

numMeas=size(zUV,2);
if(size(SR,3)==1&&numMeas>1)
    SR=repmat(SR,[1,1,numMeas]);
end

h=@(x)uv2SpherAng(x,systemType,Ms,Muv);

numCubPts=length(w);

azEl=zeros(2,numMeas);
RAzEl=zeros(2,2,numMeas);
for curMeas=1:numMeas
    xiCur=h(transformCubPoints(xi,zUV(:,curMeas),SR(:,:,curMeas)));
    meanVal=meanDirectionAzEl(xiCur,systemType,w);
    %For the computation of the covariance matrix, we need to wrap the
    %differences.
    P=zeros(2,2);
    for k=1:numCubPts
        %The azimuth is wrapped. The elevation is not wrapped.
        diffVal=xiCur(:,k)-meanVal;
        diffVal(1)=wrapRange(diffVal(1),-pi,pi);
        P=P+w(k)*(diffVal*diffVal');
    end
    azEl(:,curMeas)=meanVal;
    RAzEl(:,:,curMeas)=P;
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
