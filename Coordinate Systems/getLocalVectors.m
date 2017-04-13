function vLocal=getLocalVectors(vGlobal,uList)
%%GETGLOBALVECTORS Change a collection of local vectors into global
%                  vectors using the local coordinate axes. This multiplies
%                  the components of the vectors by the corresponding
%                  coordinate axis vectors.
%
%INPUTS:  vGlobal  A 3XN matrix of N global vectors that are to be
%                  converted.              
%           uList  A 3X3XN matrix of orthonormal unit coordinate axes that
%                  are associated with the local coordinate system in which
%                  the vectors are expressed. If all N vectors use the same
%                  local coordinate system, then a 3x3x1 matrix may be
%                  passed instead.
%
%OUTPUTS:   vLocal The vectors in the local coordinate system.
%
%January 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numVel=size(vGlobal,2);
    
    if(size(uList,3)==1)
        uList=repmat(uList,[1,1,numVel]);
    end
    
    vLocal=zeros(3,numVel);
    for curVel=1:numVel
        uTemp=uList(:,:,curVel).';
        vLocal(:,curVel)=vGlobal(1,curVel)*uTemp(:,1)+vGlobal(2,curVel)*uTemp(:,2)+vGlobal(3,curVel)*uTemp(:,3);
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
