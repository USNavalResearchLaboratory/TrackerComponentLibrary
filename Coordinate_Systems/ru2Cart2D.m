function zC=ru2Cart2D(z,useHalfRange,zTx,zRx,M)
%%RU2CART2D Convert points in 2D bistatic r-u (or r-u-v) coordinates into
%         Cartesian coordinates. If r-u coordinates are used, the target
%         is assumed to be in front of the receiver (local y coordinate is
%         positive). r-u coordinates consist of a bistatic range and
%         direction cosines at the receiver. The "direction cosine" u is
%         just the x coordinate of a unit vector from the receiver to the
%         target in the 2D coordinate system at the receiver. This
%         basically assumes that the boresight direction of the receiver is
%         the y axis. Assuming the target is in front of the
%         receiver, the second unit vector coordinate is not needed.
%         However, 2D r-u-v coordinates include the second component.
%         For monostatic coordinate conversion where the range is a one-way
%         range, set useHalfRange=true and zTx and zRx to the same value.
%
%INPUTS: z A 2XN matrix of vectors with elements [r;u], where r is the
%          bistatic range from the transmitter to the target to the
%          receiver, and u is a direction cosine. Each u value should have
%          a magnitude less than or equal to one. Alternatively, one
%          can pass a 3XN matrix of [r;u;v] vectors where [u;v] form a
%          full unit vector in the receiver's local 2D Cartesian
%          coordinates.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided, or an
%          empty matrix is passed, is false.
%      zTx The 2XN [x;y] location vectors of the transmitters in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the transmitters are assumed to be at the
%          origin. If only a single vector is passed, then the transmitter
%          location is assumed the same for all of the target states being
%          converted. zTx can have more than 2 rows; additional rows are
%          ignored.
%      zRx The 2XN [x;y] location vectors of the receivers in Cartesian
%          coordinates.  If this parameter is omitted or an empty matrix
%          is passed, then the receivers are assumed to be at the origin.
%          If only a single vector is passed, then the receiver location
%          is assumed the same for all of the target states being
%          converted. zRx can have more than 2 rows; additional rows are
%          ignored.
%        M A 2X2XN hypermatrix of the rotation matrices to go from the
%          alignment of the global coordinate system to that at the
%          receiver. The y-axis of the local coordinate system of the
%          receiver is the pointing direction of the receiver. If omitted
%          or an empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(2) --the
%          identity matrix is used. If only a single 2X2 matrix is passed,
%          then it is assumed to be the same for all of the N conversions.
%
%OUTPUTS: zC The 2XN matrix of the converted points in [x;y;z] Cartesian
%            coordinates.
%
%The 2D r-u conversion is related to the 3D r-u-v conversion, details of
%which are given in [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%February 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    N=size(z,2);

    if(nargin<5||isempty(M))
        M=repmat(eye(2),[1,1,N]);
    elseif(size(M,3)==1)
        M=repmat(M,[1,1,N]);
    end

    if(nargin<4||isempty(zRx))
        zRx=zeros(2,N);
    elseif(size(zRx,2)==1)
        zRx=repmat(zRx,[1,N]);
    end

    if(nargin<3||isempty(zTx))
        zTx=zeros(2,N);
    elseif(size(zTx,2)==1)
        zTx=repmat(zTx,[1,N]);
    end
    
    if(nargin<2||isempty(useHalfRange))
        useHalfRange=false;
    end

    %Extract the components
    rB=z(1,:);
    u=z(2,:);
    
    %The bistatic range is used in the conversions below.
    if(useHalfRange)
        rB=2*rB; 
    end
    
    %If a full unit vector is given.
    hasV=size(z,1)>2;
    
    zC=zeros(2,N);
    for curPoint=1:N
        if(hasV)
            uVec=z(2:3,curPoint);
        else
            %A unit vector pointing in the direction of the target from the
            %receiver. The real command is also if finite precision errors
            %made u slightly larger than 1.
            uVec=[u(curPoint);real(sqrt(1-u(curPoint)^2))];
        end
        %The transmitter location in the receiver's local coordinate
        %system.
        zTxL=M(:,:,curPoint)*(zTx(1:2,curPoint)-zRx(1:2,curPoint));

        r1=(rB(curPoint)^2-norm(zTxL)^2)/(2*(rB(curPoint)-dot(uVec,zTxL)));

        %This is the Cartesian location in the local coordinate system of
        %the receiver.
        zL=r1*uVec;
        %Convert to global Cartesian coordinates.
        zC(:,curPoint)=M(:,:,curPoint)\zL+zRx(1:2,curPoint);
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
