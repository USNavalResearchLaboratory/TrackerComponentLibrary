function M=randRotMat(numDim,numRotMat)
%%RANDROTMAT Generate a uniformly distributed rotation matrix in 2D or 3D.
%
%INPUTS: numDim The number of dimensions of the rotation. This can be
%               either 2 or 3.
%     numRotMat The number of rotation matrices to generate. The default if
%               omitted or an empty matrix is passed is 1.
%
%OUTPUTS: M A numDimXnumDimXnumRotMat set of rotation matrices such that
%           for each matrix M*V rotates the vector v to a uniformly
%           distributed orientation.
%
%In 2D, there is only one degree of freedom for rotations, so a single
%rotation angle can be chosen uniformely from 0->2*pi. In 3D, there are 3
%degrees of freedom to specify rotations. In n dimensions, there are
%n*(n-1)/2 degrees of freedom.
%
%In 3D, things are more difficult One cannot choose a sequence of
%uniformely distributed set of 3 Euler angles, because the resulting
%rotations will not be uniformely distributed in space. The technique
%chosen here is to choose Euler angles for a zxz rotation series as
%described in [1] where the rotations about the z-axes are uniformly
%distributed from 0-2*pi, and the rotation about the x-axis is not
%uniformly distributed.
%
%Note that general random orthonormal matrices as produced by randOrthoMat
%can have determinantes of +1 or -1, whereal all valid rotation matrices
%have determinante of +1.
%
%REFERENCES:
%[1] M. D. Shuster, "Uniform attitude probability distributions," The
%    Journal of the Astronautical Sciences, vol. 51, no. 4, pp. 451-475,
%    Oct. - Dec. 2003.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(numRotMat))
    numRotMat=1;
end

M=zeros(numDim,numDim,numRotMat);
switch(numDim)
    case 2%In 2D, there is only 1 degree of freedom.
        for k=1:numRotMat
            theta=2*pi*rand(1);
            cosT=cos(theta);
            sinT=sin(theta);
            M(:,:,k)=[cosT, -sinT;
                      sinT, cosT];
        end
    case 3
        for k=1:numRotMat
            %Sample the uniformely distributed rotations about the z-axes.
            theta1=2*pi*rand(1);
            theta3=2*pi*rand(1);
    
            %Sample the non-uniformly distributed rotation about the x-axis.
            theta2=acos(1-2*rand(1));
    
            %Turn the Euler angles into a rotation matrix.
            M(:,:,k)=Euler3Ang2RotMat(theta1,theta2,theta3,'zxz');
        end
    otherwise
        error('The dimensionality of the rotation matrix must be 2 or 3.')
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
