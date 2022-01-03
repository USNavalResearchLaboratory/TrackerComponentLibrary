function [zUV,RUV]=cameraCoords2UVCoordsCubature(zCam,SR,A,xi,w)
%%CAMERACOORDS2UVCOORDSCUBATURE Given Gaussian measurements in the local
%           coordinate system of a perspective camera, convert them to uv
%           direction cosines. This function does not then rotate them to
%           be aligned with a global set of coordinates.
%
%INPUTS: zCam A 2XnumPts set of [x;y] coordinates on the CCD of the camera
%          measuring the directions. This is in distance units, not pixels.
%       SR The 2X2XnumPts lower-triangular square roots of the covariance
%          matrices associated with zCam. If all of the matrices are the
%          same, then this can just be a single 2X2 matrix.
%        A A 3X3 matrix, as described in the comments to
%          cameraCoords2UVCoords. The third row must be [0,0,1].
%       xi A 2XnumCubaturePoints matrix of cubature points for the numeric
%          integration. If this and the final parameter are omitted or
%          empty matrices are passed, then fifthOrderCubPoints is used to
%          generate cubature points.
%        w A numCubaturePointsX1 vector of the weights associated with the
%          cubature points.
%
%OUTPUTS: zUV The 2XnumPts approximate means of the PDF of the polar
%             measurements converted to [x;y] Cartesian coordinates.
%       RCart The 2X2XnumPts set of approximate covariance matrices for the
%             numPts estimates.
%
%Details of the numerical integration used in the conversion are given in
%[1]. The perspective camera model is in [2] and more details on the model
%used are given in the comments to cameraCoords2UVCoords.
%
%EXAMPLE:
%This performs a conversion and show that the converted values are
%consistent in terms of the normalized estimation errors squared (NEES).
%That is the NEES will typically lie within the 99% confidence bouncs
%(which are also evaluated here):
% numMCRuns=1e4;
% f=35e-3;%35mm camera;
% A=diag([f,f,1]);
% SR=diag([1e-4,1e-4]);
% uvTrue=[-0.2;
%          0.5];
% zCamTrue=uvCoords2CameraCoords(uvTrue,A);
% [xi,w]=fifthOrderCubPoints(2);
% NEES=0;
% for curRun=1:numMCRuns
%     zCam=zCamTrue+SR*randn(2,1);
%     [zUV,RUV]=cameraCoords2UVCoordsCubature(zCam,SR,A,xi,w);
%     
%     diff=zUV-uvTrue;
%     NEES=NEES+diff'*inv(RUV)*diff;
% end
% NEES=NEES/(2*numMCRuns)
% NEESBounds=getNEESConfBounds(0.99,2,numMCRuns)
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%[2] J. Kannala, J. Heikkilä, and S. S. Brandt, "Geometric camera
%    calibration," in Wiley Encyclopedia of Computer Science and
%    Engineering, B. W. Wah, Ed., 2007, vol. 1.
%
%November 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numPoints=size(zCam,2);
    
    if(nargin<4||isempty(xi))
        [xi,w]=fifthOrderCubPoints(2);
    end
    
    if(size(SR,3)==1)
        SR=repmat(SR,[1,1,numPoints]);
    end
    
    h=@(z)cameraCoords2UVCoords(z,A,[],false);
    
    zUV=zeros(2,numPoints);
    RUV=zeros(2,2,numPoints);
    for curMeas=1:numPoints
        [zUV(:,curMeas), RUV(:,:,curMeas)]=calcCubPointMoments(zCam(:,curMeas),SR(:,:,curMeas),h,xi,w);
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
