function FBig=blkDiagRep(F,numTimes)
%%BLKDIAGREP Create a block diagonal matrix either where where the matrices
%            in hypermatrix F are placed on the diagonal or the matrix F is
%            repeated numTime on the main diagonal.
%
%INPUTS: F If numTimes is provided, then F is a two-dimensional matrix that
%          is repeated on the main diagonal of a larger matrix numTimes.
%          Otherwise, if the numTimes argument is omitted, F is a 3D
%          hypermatrix and each 2D element as parameterized by the third
%          will be placed on the main diagonal of a large 2D matrix.
% numTimes The number of times (>=1) that a matrix is repeated on the block
%          diagonal of the larger matrix if F is a 2D matrix to be
%          repeated. Otherwise, this parameter is ignored.
%
%OUTPUTS: FBig A large block diagonal matrix either with the 2D matrix F
%              repeated numTimes on the diagonal or a block diagonal matrix
%              with the 2D components matrices in a 3D hypermatrix F placed
%              on the main diagonal.
%
%The Matlab function blkdiag requires that the number of times a matrix is
%put onto a diagonal be fixed from the start. For example, for 2D, one uses
%blkdiag(F,F) and for 3D one uses blkdiag(F,F,F). This function allows one
%to essentially use blkdiag when the number of matrices can be variable.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

xDim=size(F,1);
yDim=size(F,2);
if(nargin==2)
    numDimX=numTimes*xDim;
    numDimY=numTimes*yDim;
    FBig=zeros(numDimX,numDimY);

    for curBlock=1:numTimes
        minIdxX=(curBlock-1)*xDim+1;
        maxIdxX=curBlock*xDim;
        spanX=minIdxX:maxIdxX;

        minIdxY=(curBlock-1)*yDim+1;
        maxIdxY=curBlock*yDim;
        spanY=minIdxY:maxIdxY;

        FBig(spanX,spanY)=F;
    end
else
    numMat=size(F,3);

    numDimX=numMat*xDim;
    numDimY=numMat*yDim;

    FBig=zeros(numDimX,numDimY);
    for curBlock=1:numMat
        minIdxX=(curBlock-1)*xDim+1;
        maxIdxX=curBlock*xDim;
        spanX=minIdxX:maxIdxX;

        minIdxY=(curBlock-1)*yDim+1;
        maxIdxY=curBlock*yDim;
        spanY=minIdxY:maxIdxY;

        FBig(spanX,spanY)=F(:,:,curBlock);
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
