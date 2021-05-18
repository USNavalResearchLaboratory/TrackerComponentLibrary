function P=clipCovMat(P,minVals,maxVals,isSqrt)
%%CLIPCOVMAT Force the diagonal elements of a covariance matrix to be no
%         smaller than and no larger than certain values while still
%         keeping the matrix a valid covariance matrix. This function can
%         also be used with lower-triangular square root covariance
%         matrices. This function is useful if one wishes to set a maximum
%         value on the uncertainty in certain components of a target state
%         based on reasonable physical lower and upper bounds. For example,
%         if a target type cannot go over 100m/s, but the 99% confidence
%         region of the target covers +/-10^8 meters per second, it is
%         reasonable to shrink the region.
%
%INPUTS: P An NXNXnumMats set of numMats covariance matrices (or lower-
%          triangular square root covariance matrices if isSqrt is true).
%          If minVals is passed, then it is assumed that none of the
%          diagonal elements equals zero.
%  minVals A 1XN or NX1 vector that contains the minimum positive value
%          allowed on each diagonal element of the matrices. If only a
%          scalar is passed, then the same value is used for all N
%          dimensions. If no lower bounds are set, then an empty matrix can
%          be passed.
%  maxVals A 1XN or NX1 vector that contains the maximum positive value
%          allowed on each diagonal element of the matrices. If only a
%          scalar is passed, then the same value is used for all N
%          dimensions. If no upper bounds are set, then an empty matrix can
%          be passed. it is assumed that maxVals(i)>minVals(i).
%   isSqrt A boolean value indicating whether P holds square root
%          covariance matrices or lower-triangular square root matrices (as
%          one might get with chol using the 'lower' option). The default
%          if omitted or an empty matrix is passed is false.
%
%OUTPUTS: P The NXNXnumMats set of matrices after scaling so the diagonal
%           elements are no larger than the values in maxVals and no
%           smaller than the values in minVals.
%
%The diagonal elements of a covariance matrix are related to the largest
%rectangle that a certain uncertainty ellipsoid can be placed in. This is
%derived in Appendix A of [1].
%
%We cannot just shrink the diagonal elements without adjusting the cross
%terms or the matrix can fail to remain positive definite. It is well known
%that multiplying a random variable v by a matrix S, S*v, changes the
%covariance matrix of the random variable from P to S*P*S'. Thus, we choose
%a diagonal matrix S such that the diagonal elements of P that are too
%large are scaled to the appropriate size. One would hope that this type of
%scaling would retain reasonable covariance cross terms for a target
%tracking algorithm. Due to the scaling, this function cannot handle zero
%diagonal elements that are too small.
%
%REFERENCES:
%[1] A. B. Poore, "Complexity reduction in MHT/MFA tracking," in
%    Proceedings of SPIE: Signal and Data Processing of Small Targets, vol.
%    5913, San Diego, CA, 31 Jul. 2005.
%
%March 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numDim=size(P,1);
numMats=size(P,3);

if(nargin<4||isempty(isSqrt))
    isSqrt=false;
end

if(nargin<3)
    maxVals=[];%isscalar will be false. 
end

if(isscalar(minVals))
    minVals=minVals*ones(numDim,1);
end

if(isscalar(maxVals))
    maxVals=maxVals*ones(numDim,1);
end

if(isempty(maxVals)&&isempty(minVals))
    return; 
end

if(isSqrt)
    %If P is actually lower-triangular square root covariance matrices
    %rather than the matrices themselves. 
    for curMat=1:numMats
        PDiagCur=diag(P(:,:,curMat)*P(:,:,curMat)');
        
        S=eye(numDim,numDim);
        if(~isempty(maxVals))
            for curDim=1:numDim
               if(PDiagCur(curDim)>maxVals(curDim))
                   S(curDim,curDim)=sqrt(maxVals(curDim)/PDiagCur(curDim));
               end
            end
        end
        if(~isempty(minVals))
            for curDim=1:numDim
               if(PDiagCur(curDim)<minVals(curDim))
                   S(curDim,curDim)=sqrt(minVals(curDim)/PDiagCur(curDim));
               end
            end
        end
        P(:,:,curMat)=S*P(:,:,curMat);
    end
else
    for curMat=1:numMats
        S=eye(numDim,numDim);
        if(~isempty(maxVals))
            for curDim=1:numDim
               if(P(curDim,curDim,curMat)>maxVals(curDim))
                   S(curDim,curDim)=sqrt(maxVals(curDim)/P(curDim,curDim,curMat));
               end
            end
        end
        if(~isempty(minVals))
            for curDim=1:numDim
               if(P(curDim,curDim,curMat)<minVals(curDim))
                   S(curDim,curDim)=sqrt(minVals(curDim)/P(curDim,curDim,curMat));
               end
            end
        end
        P(:,:,curMat)=S*P(:,:,curMat)*S';
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
