function u=getUDirection2D(zC,zRx,M,includeV)
%%GETUDIRECTION2D Convert points in 2D Cartesian coordinates into a local
%         direction cosine component u in 2D. The "direction cosine" u is
%         just the x coordinate of a unit vector from the receiver to the
%         target in the coordinate system at the receiver. This basically
%         assumes that the boresight direction of the receiver is the y
%         axis. Assuming the target is in front of the receiver, the second
%         unit vector coordinate is not needed. However, with the includeV
%         option, it can be returned.
%
%INPUT: zC A 2XN matrix of Cartesian points in global [x;y] Cartesian
%          coordinates. Extra rows are ignored.
%      zRx The 2XN [x;y] location vectors of the receiver in Cartesian
%          coordinates.  If this parameter is omitted or an empty matrix
%          is passed, then the receivers are assumed to be at the origin.
%          If only a single vector is passed, then the receiver location
%          is assumed the same for all of the target states being converted
%          zRx can have more than 2 rows; additional rows are ignored.
%        M A 2X2XN hypermatrix of the rotation matrices to go from the
%          alignment of the global coordinate system to that at the
%          receiver. The y-axis of the local coordinate system of the
%          receiver is the pointing direction of the receiver. If omitted
%          or an empty matrix is passed, then it is assumed that the local
%          coordinate system is aligned with the global and M=eye(2) --the
%          identity matrix is used. If only a single 2X2 matrix is passed,
%          then it is assumed to be the same for all of the N conversions.
% includeV An optional boolean value indicating whether a second direction
%          cosine component should be included. The u direction cosine is
%          one parts of a 2D unit vector. Generally, one might assume that
%          the target is in front of the sensor, so the second component
%          would be positive and is not needed. However, the second
%          component can be included if ambiguity exists. The default if
%          this parameter is omitted or an empty matrix is passed is false.
%
%OUTPUT: u The 1XN (or 2XN if includeV is true) matrix of direction cosine
%          values u (or [u;v]) for each of the points in zC.
%
%The conversion is related to the 3D r-u-v conversion, details of which are
%given in [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%June 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release

    if(nargin<4||isempty(includeV))
        includeV=false;
    end
 
    N=size(zC,2);
    
    if(nargin<3||isempty(M))
        M=repmat(eye(2),[1,1,N]);
    elseif(size(M,3)==1)
        M=repmat(M,[1,1,N]);
    end

    if(nargin<2||isempty(zRx))
        zRx=zeros(2,N);
    elseif(size(zRx,2)==1)
        zRx=repmat(zRx,[1,N]);
    end

    if(includeV==true)
        u=zeros(2,N);
    else
        u=zeros(1,N);
    end
    
    for curPoint=1:N
        %The target location in the receiver's coordinate system.
        zCL=M(:,:,curPoint)*(zC(1:2,curPoint)-zRx(1:2,curPoint));
        r1=norm(zCL);%Receiver to target.
        u(1,curPoint)=zCL(1)/r1;
        if(includeV)
           u(2,curPoint)=zCL(2)/r1; 
        end
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
