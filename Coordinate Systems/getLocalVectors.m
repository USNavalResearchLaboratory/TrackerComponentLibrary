function vLocal=getLocalVectors(vGlobal,uList)
%%GETLOCALVECTORS Change a collection of global vectors into local
%                 vectors using the global coordinate axes. This multiplies
%                 the components of the vectors by the corresponding
%                 coordinate axis vectors. Given the global coordinate axes
%                 as orthonormal unit vectors, the inverse basis vectors
%                 are obtained as the transposes of each 3X3 matrix in
%                 uList.
%
%INPUTS: vGlobal A numDimsXN matrix of N global vectors that are to be
%                converted.              
%          uList A numDimsXnumDimsXN matrix of orthonormal unit coordinate
%                axes that are associated with the local coordinate system
%                in which the vectors are expressed. If all N vectors use
%                the same local coordinate system, then a numDimsxnumDimsx1
%                matrix can be passed instead.
%
%OUTPUTS: vLocal The vectors in the local coordinate system.
%
%January 2015 David Karnick, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    numVel=size(vGlobal,2);
    
    if(size(uList,3)==1)
        uList=repmat(uList,[1,1,numVel]);
    end
    
    numDims=size(vGlobal,1);
    vLocal=zeros(numDims,numVel);
    for curVel=1:numVel
        uTemp=uList(:,:,curVel).';
        for curDim=1:numDims
            vLocal(:,curVel)=vLocal(:,curVel)+vGlobal(curDim,curVel)*uTemp(:,curDim);
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
