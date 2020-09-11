function zC=ruv2Cart(z,useHalfRange,zTx,zRx,M)
%RUV2CART Convert points in bistatic r-u-v (or r-u-v-w) coordinates into
%         Cartesian coordinates. If r-u-v coordinates are used, the target
%         is assumed to be in front of the receiver (local z coordinate is
%         positive). r-u-v coordinates consist of a bistatic range and
%         direction cosines at the receiver. The "direction cosines" u and
%         v are just the x and y coordinates of a unit vector from the
%         receiver to the target in the coordinate system at the receiver.
%         This basically assumes that the boresight direction of the
%         receiver is the z axis. Assuming the target is in front of the
%         receiver, the third unit vector coordinate is not needed.
%         However, r-u-v-w coordinates include the third component.
%         For monostatic coordinate conversion where the range is a one-way
%         range, set useHalfRange=true and zTx and zRx to the same value.
%
%INPUTS: z A 3XN matrix of vectors with elements [r;u;v], where r is the
%          bistatic range from the transmitter to the target to the
%          receiver, and u and v are direction cosines. Each u,v pair
%          should have a magnitude less than or equal to one. If the
%          magnitude is greater than one, then the pair is normalized
%          before conversion to avoid imaginary results. Alternatively, one
%          can pass a 4XN matrix of [r;u;v;w] vectors where [u;v;w] form a
%          full unit vector in the receiver's local 3D Cartesian
%          coordinates.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided, or an
%          empty matrix is passed, is false.
%      zTx The 3XN [x;y;z] location vectors of the transmitters in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the transmitters are assumed to be at the
%          origin. If only a single vector is passed, then the transmitter
%          location is assumed the same for all of the target states being
%          converted. zTx can have more than 3 rows; additional rows are
%          ignored.
%      zRx The 3XN [x;y;z] location vectors of the receivers in Cartesian
%          coordinates.  If this parameter is omitted or an empty matrix
%          is passed, then the receivers are assumed to be at the origin.
%          If only a single vector is passed, then the receiver location
%          is assumed the same for all of the target states being
%          converted. zRx can have more than 3 rows; additional rows are
%          ignored.
%        M A 3X3XN hypermatrix of the rotation matrices to go from the
%          alignment of the global coordinate system to that at the
%          receiver. The z-axis of the local coordinate system of the
%          receiver is the pointing direction of the receiver. If omitted
%          or an empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(3) --the
%          identity matrix is used. If only a single 3X3 matrix is passed,
%          then it is assumed to be the same for all of the N conversions.
%
%OUTPUTS: zC The 3XN matrix of the converted points in [x;y;z] Cartesian
%            coordinates.
%
%Basic u and v direction cosines do not specify which side of the radar the
%target is on. That is, they do not specify the sign of z. This
%performs the conversion assuming that z is positive in the local
%coordinate system of the receiver. Details of the conversion are given in
%[1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    N=size(z,2);

    if(nargin<5||isempty(M))
        M=repmat(eye(3),[1,1,N]);
    elseif(size(M,3)==1)
        M=repmat(M,[1,1,N]);
    end

    if(nargin<4||isempty(zRx))
        zRx=zeros(3,N);
    elseif(size(zRx,2)==1)
        zRx=repmat(zRx,[1,N]);
    end

    if(nargin<3||isempty(zTx))
        zTx=zeros(3,N);
    elseif(size(zTx,2)==1)
        zTx=repmat(zTx,[1,N]);
    end
    
    if(nargin<2||isempty(useHalfRange))
        useHalfRange=false;
    end

    %Extract the components
    rB=z(1,:);
    u=z(2,:);
    v=z(3,:);
    
    %The bistatic range is used in the conversions below.
    if(useHalfRange)
        rB=2*rB; 
    end
    
    %If a full unit vector is given.
    hasW=size(z,1)>3;
    
    zC=zeros(3,N);
    for curPoint=1:N
        if(hasW)
            uVec=z(2:4,curPoint);
        else
            %If the magnitude is too large, normalize it so that one does
            %not encounter imaginary numbers.
            uvMag2=u(curPoint)^2+v(curPoint)^2;
            if(uvMag2>1)
                uvMag=sqrt(uvMag2);
                u(curPoint)=u(curPoint)/uvMag;
                v(curPoint)=v(curPoint)/uvMag;
            end

            %A unit vector pointing in the direction of the target from the
            %receiver. The real command is also if the normalization didn't
            %quite work due to finite precision errors.
            uVec=[u(curPoint);v(curPoint);real(sqrt(1-u(curPoint)^2-v(curPoint)^2))];
        end
        %The transmitter location in the receiver's local coordinate
        %system.
        zTxL=M(:,:,curPoint)*(zTx(1:3,curPoint)-zRx(1:3,curPoint));

        r1=(rB(curPoint)^2-norm(zTxL)^2)/(2*(rB(curPoint)-dot(uVec,zTxL)));
        
        %This deals with the case where a point with zero range in the
        %monostatic case is passed. In the bistatic case, a zero range is
        %invalid.
        if(rB(curPoint)==0)
            r1=0;
        end

        %This is the Cartesian location in the local coordinate system of
        %the receiver.
        zL=r1*uVec;
        %Convert to global Cartesian coordinates.
        zC(:,curPoint)=M(:,:,curPoint)\zL+zRx(1:3,curPoint);
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
