function z=getUDirection3D(zC,zRx,M)
%%GETUDIRECTION3D Convert a position into a single local direction cosines
%         u. Direction cosines u and v are just the x and 
%         coordinates of a unit vector from the receiver to the target in
%         the coordinate system at the receiver. This assumes that the
%         boresight direction of the receiver is the local z axis.
%
%INPUTS: zC A 3XN matrix of Cartesian points in global [x;y;z] Cartesian
%          coordinates.
%      zRx The 3X1 [x;y;z] location vector of the receiver in Cartesian
%          coordinates.  If this parameter is omitted or an empty matrix is
%          passed, then the receiver is assumed to be at the origin.
%        M A 3X3 rotation matrix to go from the alignment of the global
%          coordinate system to that at the receiver. The z-axis of the
%          local coordinate system of the receiver is the pointing
%          direction of the receiver. If omitted or an empty matrix is
%          passed, then it is assumed that the local coordinate system is
%          aligned with the global and M=eye(3) --the identity matrix is
%          used.
%
%OUTPUTS: z The 1XN vector of the v direction cosines for the points in zC.
%
%Details of the conversion are given in [1].
%
%REFERENCES:
%[1] David F. Crouse , "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%March 2021 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

N=size(zC,2);

if(nargin<3||isempty(M))
    M=eye(3);
end

if(nargin<2||isempty(zRx))
    zRx=zeros(3,1);
end

%Allocate space for the return values.
z=zeros(1,N);
for curPoint=1:N
    %The target location in the receiver's coordinate system.
    zCL=M*(zC(:,curPoint)-zRx);

    %Perform the conversion.
    r1=norm(zCL);%Receiver to target.
    z(1,curPoint)=zCL(1)/r1;
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
