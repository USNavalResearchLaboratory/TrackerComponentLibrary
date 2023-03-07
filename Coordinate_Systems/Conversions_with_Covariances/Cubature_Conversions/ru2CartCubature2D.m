function [zCart,RCart]=ru2CartCubature2D(z,SR,useHalfRange,zTx,zRx,M,xi,w)
%%RU2CARTCUBATURE2D Use cubature integration to approximate the moments of
%         a  noise-corrupted bistatic r-u (or r-u-w) location into 2D
%         Cartesian coordinates. If r-u coordinates are used, the target
%         is assumed to be in front of the receiver (local y coordinate is
%         positive). r-u coordinates consist of a bistatic range and a
%         direction cosine at the receiver. The "direction cosine" u is
%         just the x coordinate of a unit vector from the receiver to the
%         target in the coordinate system at the receiver. This basically
%         assumes that the boresight direction of the receiver is the y
%         axis. Assuming the target is in front of the receiver, the second
%         unit vector coordinate is not needed. However, r-u-w coordinates
%         include the third component. For monostatic coordinate conversion
%         where the range is a one-way range, set useHalfRange=true and zTx
%         and zRx to the same value.
%
%INPUTS: z A 2XN matrix of vectors with elements [r;u], where r is the
%          bistatic range from the transmitter to the target to the
%          receiver, and u is the direction cosine, which should have a
%          magnitude less than or equal to one. If the magnitude is greater
%          than one, then it is normalized before conversion to avoid
%          imaginary results.
%       SR The 2X2XN lower-triangular square roots of the covariance
%          matrices associated with z. If all of the matrices are
%          the same, then this can just be a single 2X2 matrix.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided, or an
%          empty matrix is passed, is false.
%      zTx The 2X1 [x;y] location vector of the transmitter in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the transmitter is assumed to be at the
%          origin. zTx can have more than 3 rows; additional rows are
%          ignored.
%      zRx The 2X1 [x;y] location vector of the receiver in Cartesian
%          coordinates.  If this parameter is omitted or an empty matrix is
%          passed, then the receiver is assumed to be at the origin. zRx
%          can have more than 3 rows; additional rows are ignored.
%        M A 2X2 rotation matrices to go from the alignment of the global
%          coordinate system to that at the receiver. If omitted or an
%          empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(2) --the
%          identity matrix is used.
%       xi A 2XnumCubaturePoints (for position-only) matrix of cubature
%          points for the numeric integration. If this and the final
%          parameter are omitted or empty matrices are passed, then
%          fifthOrderCubPoints is used to generate cubature points.
%        w A numCubaturePointsX1 vector of the weights associated with the
%          cubature points.
%
%OUTPUTS: zCart The 2XN approximate means of the PDF of the polar
%               measurements converted to [x;y] Cartesian coordinates.
%         RCart The 2X2XN set of approximate covariance matrices for the N
%               estimates.
%
%Details of the numerical integration used in the conversion are given in
%[1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%September 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numPoints=size(z,2);
    
    if(nargin<7||isempty(xi))
        [xi,w]=fifthOrderCubPoints(2);
    end

    if(nargin<6||isempty(M))
        M=eye(2);
    end

    if(nargin<5||isempty(zRx))
        zRx=zeros(2,1);
    end

    if(nargin<4||isempty(zTx))
        zTx=zeros(2,1);
    end

    if(nargin<3||isempty(useHalfRange))
        useHalfRange=false;
    end
    
    numMeas=size(z,2);

    if(size(SR,3)==1)
        SR=repmat(SR,[1,1,numMeas]);
    end
    
    h=@(z)ru2Cart2D(z,useHalfRange,zTx(1:2,:),zRx(1:2,:),M);
    
    zCart=zeros(2,numPoints);
    RCart=zeros(2,2,numPoints);
    for curMeas=1:numPoints
        [zCart(:,curMeas), RCart(:,:,curMeas)]=calcCubPointMoments(z(:,curMeas),SR(:,:,curMeas),h,xi,w);
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
