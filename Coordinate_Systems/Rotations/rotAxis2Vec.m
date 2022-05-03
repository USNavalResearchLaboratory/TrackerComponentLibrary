function R=rotAxis2Vec(u,axisVal,method)
%%ROTAXIS2VEC Get a rotation matrix that rotates the a coordinate axis
%             (e.g. the z-axis) in the direction of a given vector u. That
%             is, u=R*z. Another way to look at it for non-unit vectors
%             when considering three dimensions is R'*u=[0;0;norm(u)]
%             The solution to the problem is not unique, so two methods
%             are given. This algorithm works with the number of
%             dimensions>=1.
%
%INPUTS: u1 A numDimXN unit vector representing the direction into which
%           the chosen axis is to be rotated. Non-unit vectors that are
%           passed will be normalized. For N>1, multiple rotation matrices
%           are returned.
%   axisVal The axis that is chosen. This selection the index can be a
%           number from 1 to numDim, or it can be 'x', 'y' or 'z' to select
%           dimensions 1, 2, and 3. if this parameter is omitted and
%           numDim=3, then axis=3 is used, otherwise axis=1 is used.
%    method This selects the type of rotation used. Possible values are:
%           0 (The default if omitted or an empty matrix is passed). Use a
%             Householder rotation via the function HouseholderVec. This
%             has the property that R*R is the identity matrix. This method
%             is also numerically stabler than method 1.
%           1 This method only works for 3D vectors. Use the approach of
%             Appendix C of [1] to obtain the rotation matrix corresponding
%             to the shortest rotation angle between the axis and the
%             point. This method has a loss of precision when the axis is
%             very close to u.
%
%OUTPUTS: R A numDumXnumDimXN set of rotation matrices such that
%           u=R*theAxis (within a scaling factor if u was not a unit
%           vector) for all of the vectors in u, where theAxis is the
%           chosen axis. That is [1;0;0] for 'x' (or 1), [0;1;0] for 'y' or
%           2, and [0;0;1] for z or 3.
%
%EXAMPLE:
% uVec=[53;183;-225;86;31;-130;-43;34];
% axisVal=5;
% R=rotAxis2Vec(uVec,axisVal);
% R'*uVec
%One will see that the rotated uVec is about 
%[0;0;0;0;339.3891571632777;0;0;0];
%
%REFERENCES:
%[1] D. F. Crouse, "On measurement-based light-time corrections for
%    bistatic orbital debris tracking," IEEE Transactions on Aerospace and
%    Electronics Systems, vol. 51, no. 3, pp. 2502-2518,  Jul. 2015.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
    numDim=size(u,1);
    numVec=size(u,2);
    
    if(nargin<3||isempty(method))
        method=0; 
    end
    
    if(method==1&&numDim~=3)
        error('This method only works for 3D vectors.') 
    end

    if(nargin<2||isempty(axisVal))
        if(numDim==3)
            axisVal=3;
        else
            axisVal=1;
        end
    end
    
    if(ischar(axisVal))
        switch(axisVal)
            case 'x'
                axisVal=1;
            case 'y'
                axisVal=2;
            case 'z'
                axisVal=3;
            otherwise
                error('Invalid axis specified.')
        end
    end
    
    if(axisVal~=fix(axisVal)||~isreal(axisVal)||axisVal<1||axisVal>numDim)
        error('Invalid axis specified.')
    end

    %The scalar special case.
    if(numDim==1)
        R=ones(1,1,numVec);
        return;
    end

    R=zeros(numDim,numDim,numVec);
    switch(method)
        case 0%Use a Householder transformation.
            %The Householder transformation puts everything in the first
            %index. The permutation vector is used to make it work for any
            %axis value.
            permVec=1:numDim;
            permVec(axisVal)=1;
            permVec(1)=axisVal;

            for curVec=1:numVec
                %The force sign option is true to make it be positive.
                [v,beta]=HouseholderVec(u(permVec,curVec),true);

                R(:,:,curVec)=eye(numDim,numDim)-beta*(v*v');
                R(:,:,curVec)=R(permVec,permVec,curVec);
            end
        case 1
            %Use the shortest rotation (not accurate for very small
            %rotations).
 
            %Normalize the vectors
            u=bsxfun(@rdivide,u,sqrt(sum(u.*u,1)));
            switch(axisVal)
                case 1
                    for curVec=1:numVec
                        %The calls to max are to deal with possible finite
                        %precision issues.
                        s=-sqrt(max(0,(1/2)*(1-u(1,curVec))));
                        dist=sqrt(u(2,curVec)^2+u(3,curVec)^2);

                        q=[sqrt(max(0,1-s^2)),0,-s*u(3,curVec)/dist,s*u(2,curVec)/dist];

                        %Deal with the case where the x axis is passed.
                        q(~isfinite(q))=0;

                        R(:,:,curVec)=quat2RotMat(q,'left');
                    end
                case 2
                    for curVec=1:numVec
                        %The calls to max are to deal with possible finite
                        %precision issues.
                        s=-sqrt(max(0,(1/2)*(1-u(2,curVec))));
                        dist=sqrt(u(1,curVec)^2+u(3,curVec)^2);
                        q=[sqrt(max(0,1-s^2)),s*u(3,curVec)/dist,0,-s*u(1,curVec)/dist];

                        %Deal with the case where the y axis is passed.
                        q(~isfinite(q))=0;

                        R(:,:,curVec)=quat2RotMat(q,'left');
                    end
                case 3
                    for curVec=1:numVec
                        %The calls to max are to deal with possible finite
                        %precision issues.
                        s=-sqrt(max(0,(1/2)*(1-u(3,curVec))));
                        dist=sqrt(u(1,curVec)^2+u(2,curVec)^2);
                        q=[sqrt(max(0,1-s^2)),-s*u(2,curVec)/dist,s*u(1,curVec)/dist,0];

                        %Deal with the case where the z axis is passed.
                        q(~isfinite(q))=0;

                        R(:,:,curVec)=quat2RotMat(q,'left');
                    end
                otherwise
                    error('Unknown axis specified')
            end
        otherwise
            error('Invalid method specified.')
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
