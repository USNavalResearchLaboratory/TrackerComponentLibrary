function R=rotAxis2Vec(u,axis)
%%ROTAXIS2VEC  Get a rotation matrix that rotates the a coordinate axis
%              (e.g. the z-axis) in the direction of a given vector u. That
%              is, u=R*z. As the solution to the problem is not unique, the
%              shortest-distance rotation is chosen.
%
%INPUTS: u1 A 3XN unit vector representing the direction into which the
%           chosen axis is to be rotated. Non-unit vectors that are passed
%           will be normalized. For N>1, multiple rotation matrices are
%           returned.
%      axis The axis that is chosen. This can be 'x', 'y', or 'z'.
%
%OUTPUTS: R A 3X3XN set of rotation matrices such that u=R*theAxis (within
%           a scaling factor if u was not a unit vector) for all of the
%           vectors in u, where theAxis is the chosen axis. That is [1;0;0]
%           for x, [0;1;0] for y, and [0;0;1] for z.
%
%In Appendix C of [1], the minimal rotation to point the z-axis in a
%desired direction is derived. A similar derivation leads to the rotations
%for the x and y axes.
%
%REFERENCES:
%[1] D. F. Crouse, "On measurement-based light-time corrections for
%    bistatic orbital debris tracking," IEEE Transactions on Aerospace and
%    Electronics Systems, vol. 51, no. 3, pp. 2502-2518,  Jul. 2015.
%
%August 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
    
    if(nargin<2||isempty(axis))
       axis='z'; 
    end

%Normalize the vectors
    u=bsxfun(@rdivide,u,sqrt(sum(u.*u,1)));
    
    numVec=size(u,2);
    R=zeros(3,3,numVec);

    switch(axis)
        case 'x'
            for curVec=1:numVec
                s=-sqrt(0.5*(1-u(1,curVec)));
                dist=sqrt(u(2,curVec)^2+u(3,curVec)^2);
                
                q=[sqrt(1-s^2),0,-s*u(3,curVec)/dist,s*u(2,curVec)/dist];
                
                %Deal with the case where the x axis is passed.
                q(~isfinite(q))=0;

                R(:,:,curVec)=quat2RotMat(q,'left');
            end
        case 'y'
            for curVec=1:numVec
                s=-sqrt(0.5*(1-u(2,curVec)));
                dist=sqrt(u(1,curVec)^2+u(3,curVec)^2);
                q=[sqrt(1-s^2),s*u(3,curVec)/dist,0,-s*u(1,curVec)/dist];
                
                %Deal with the case where the y axis is passed.
                q(~isfinite(q))=0;
                
                R(:,:,curVec)=quat2RotMat(q,'left');
            end
        case 'z'
            for curVec=1:numVec
                s=-sqrt(0.5*(1-u(3,curVec)));
                dist=sqrt(u(1,curVec)^2+u(2,curVec)^2);
                q=[sqrt(1-s^2),-s*u(2,curVec)/dist,s*u(1,curVec)/dist,0];
                
                %Deal with the case where the z axis is passed.
                q(~isfinite(q))=0;
                
                R(:,:,curVec)=quat2RotMat(q,'left');
            end
        otherwise
            error('Unknown axis specified')
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
