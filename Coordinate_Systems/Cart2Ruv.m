function z=Cart2Ruv(zC,useHalfRange,zTx,zRx,M,includeW)
%CART2RUV Convert points in Cartesian coordinates (either the global
%         system or a local system at the receiver) into local bistatic
%         r-u-v coordinates of the receiver. r-u-v coordinates consist of a
%         bistatic range and direction cosines at the receiver. The
%         "direction cosines" u and v are just the x and y coordinates of a
%         unit vector from the receiver to the target in the coordinate
%         system at the receiver. This basically assumes that the boresight
%         direction of the receiver is the z axis. Assuming the target is
%         in front of the receiver, the third unit vector coordinate is not
%         needed. However, with the includeW option, it can be provided,
%         resulting in r-u-v-w coordinates.
%
%INPUT: zC A 3XN matrix of Cartesian points in global [x;y;z] Cartesian
%          coordinates. Extra rows are ignored.
% useHalfRange A boolean value specifying whether the bistatic range value
%          should be divided by two. This normally comes up when operating
%          in monostatic mode, so that the range reported is a one-way
%          range. The default if this parameter is not provided (or an
%          empty matrix is provided) is false.
%      zTx The 3XN [x;y;z] location vectors of the transmitters in global
%          Cartesian coordinates. If this parameter is omitted or an empty
%          matrix is passed, then the transmitters are assumed to be at the
%          origin. If only a single vector is passed, then the transmitter
%          location is assumed the same for all of the points being
%          converted. zTx can have more than 3 rows; additional rows are
%          ignored.
%      zRx The 3XN [x;y;z] location vectors of the receivers in Cartesian
%          coordinates.  If this parameter is omitted or an empty matrix
%          is passed, then the receivers are assumed to be at the origin.
%          If only a single vector is passed, then the receiver location
%          is assumed the same for all of the points being converted. zRx
%          can have more than 3 rows; additional rows are ignored.
%        M A 3X3XN hypermatrix of the rotation matrices to go from the
%          alignment of the global coordinate system to that at the
%          receiver. The z-axis of the local coordinate system of the
%          receiver is the pointing direction of the receiver. If omitted
%          or an empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(3) --the
%          identity matrix is used. If only a single 3X3 matrix is passed,
%          then it is assumed to be the same for all of the N conversions.
% includeW An optional boolean value indicating whether a third direction
%          cosine component should be included. The u and v direction
%          cosines are two parts of a 3D unit vector. Generally, one might
%          assume that the target is in front of the sensor, so the third
%          component would be positive and is not needed. However, the
%          third component can be included if ambiguity exists. The default
%          if this parameter is omitted or an empty matrix is passed is
%          false.
%
%OUTPUT: z The 3XN (or 4XN if includeW is true) matrix of location vectors
%          of the points in bistatic [r;u;v] coordinates. If
%          useHalfRange=true, then the r component is half the bistatic
%          range (half the round-trip range for a monostatic scenario).
% 
%Details of the conversion are given in [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<6||isempty(includeW))
        includeW=false;
    end

    N=size(zC,2);

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
    
    %Allocate space for the return values.
    if(includeW==true)
        z=zeros(4,N);
    else
        z=zeros(3,N);
    end
    
    for curPoint=1:N
        %The target location in the receiver's coordinate system.
        zCL=M(:,:,curPoint)*(zC(1:3,curPoint)-zRx(1:3,curPoint));
        %The transmitter location in the receiver's local coordinate
        %system.
        zTxL=M(:,:,curPoint)*(zTx(1:3,curPoint)-zRx(1:3,curPoint));

    %Perform the conversion.
        r1=norm(zCL);%Receiver to target.
        r2=norm(zCL-zTxL);%Target to transmitter.

        r=r1+r2;

        u=zCL(1)/r1;
        v=zCL(2)/r1;

        z(1:3,curPoint)=[r;u;v];
        if(includeW)
            z(4,curPoint)=zCL(3)/r1;
        end
    end
    
    if(useHalfRange)
        z(1,:)=z(1,:)/2;
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
