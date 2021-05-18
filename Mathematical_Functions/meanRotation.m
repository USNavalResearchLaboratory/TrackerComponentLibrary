function meanRot=meanRotation(param1,w)
%%MEANROTATION Perform a weighted or unweighted average of a set of
%              rotation matrices or unit quaternions (which represent
%              directions). The weights are scalar and the cost function
%              finds the rotation matrix REst to minimize the sum of 
%              w(i)*norm(REst-R{i},'fro'), where R{i} is the ith rotation
%              matrix to be averaged (either given explicitly as a
%              rotation matrix or given as a quaternion).
%
%INPUTS: param1 A 3X3XN hypermatrix of N rotation matrices that are to be
%               averaged or a 4XN set of unit-quaternions that are to be
%               averaged. The function automatically detects whether
%               rotation matrices or quaternions were passed. The
%               handedness of the quaternions does not matter.
%             w An optional set of N weights for averaging the rotations.
%               If omitted, all of the rotation matrices/ unit quaternions
%               in param1 are equally weighted.
%
%OUTPUTS: meanRot The weighted mean 3X3 rotation matrix or 4X1 rotation
%                 quaternion, depending on the data type passed in param1.
%
%The algorithm of [1] for performing a weighted average of quaternions or
%rotation matrices is used. Note that the paper contains an additional
%technique for averaging quaternions with non-scalar weights, which is not
%implemented here.
%
%A quaternion can represent an orientation (rotation) in space as a 4X1
%unit vector. However, one cannot simply take the weighted average of a
%set of quaternions to get an average quaternion (that is then normalized)
%to represent the average rotation, because the quaternion representation
%is not unique. Thus, while q and -q represent the same rotation, they
%affect such an average differently, which is bad. This is why the more
%complicated algorithms of the above paper are used.
%
%REFERENCES:
%[1] F. L. Markley, Y. Cheng, J. L. Crassidis, and Y. Oshman, "Averaging
%    quaternions," Journal of Guidance, Control, and Dynamics, vol. 30,
%    no. 4, pp. 1193-1196, Jul. - Aug. 2007.
%
%September 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(size(param1,1)==3)
    usingQuaternions=false;
    RotMats=param1;
    
    numVals=size(param1,3);
else
    usingQuaternions=true;
    q=param1;
    
    numVals=size(param1,2);
end

if(nargin<2)
   w=ones(numVals,1); 
end

%If quaternions are being averaged rather than rotation matrices.
if(usingQuaternions)
    %Build M as in Equation 2 of [1].
    M=zeros(4,4);
    for curQuat=1:numVals
        M=M+w(curQuat)*(q(:,curQuat)*q(:,curQuat)');
    end

    %The solution to the cost function in Equation 13 is the eigenvector
    %corresponding to the largest eigenvalue in M.
    [V,D]=eig(M);
    %Get the indices of the sorted eigenvalues (and thus of the
    %corresponding eigenvectors).
    [~,idx]=sort(diag(D),'descend');
    %The rotation quaternion is the eigenvector associated with the
    %largest eigenvalue.
    meanRot=V(:,idx(1));%A quaternion
else%If rotation matrices are being averaged rather than quaternions.
    %Build B as in Equation 5 in [1].
    B=zeros(3,3);
    for curRot=1:numVals
       B=B+w(curRot)*RotMats(:,:,curRot); 
    end

    %This solves the equivalent formulation of Wahba's problem for
    %Equation 4.
    meanRot=findTransParam(B,eye(3));%A rotation matrix.
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
