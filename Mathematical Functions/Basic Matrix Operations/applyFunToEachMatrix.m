function X=applyFunToEachMatrix(f,A,B,matDims)
%%APPLYFUNTOEACHMATRIX Given a hypermatrix containing numerous 2D matrices,
%         apply a function to each matrix. Alternatively, the function can
%         take two parameters and be applied to corresponding matrices in
%         two sets. This function is useful, for example, if one wishes to
%         invert a collection of matrices that are given in a hypermatrix.
%         Such an operation cannot be done with something like bsxfun,
%         which operates element-by-element.
%
%INPUTS: f A function handle that takes a matrix and returns a matrix of
%          the same size (or matDims in size), or a function that takes two
%          matrices and returns a matrix, if B is provided to this
%          function.
%        A A aDim1XaDim2X... hypermatrix. Each individual matrix that
%          should be passed to f is xDim1XaDim2 in size.
%        B An optional bDim1XbDim2X... hypermatrix. There must be the same
%          number of matrices in this hypermatrix as in A, though they need
%          not be the same size. If this matrix is provided and is not
%          empty, then f will be called with two arguments. Otherwise, f
%          will be called with only one argument.
%  matDims This is a 1X2 vector that specifies the dimensions of the output
%          matrix of f. If this parameter is omitted or an empty matrix is
%          passed, then it will be assumed that the output of f has the
%          same number of rows and columns as A.
%
%OUTPUTS: X A hypermatrix with size of the rows and columns as specified
%           above and the size of the higher dimensions equal to those of
%           A. These matrices are from the output of f.
%
%For GPUs, a similar function named pagefun exists. However, a CPU-only
%counterpart does not currently exist (2017), hence the need for this
%function.
%
%EXAMPLE 1:
%Here, we have a number of matrices in a hypermatrix, and we want to invert
%all of them.
% f=@(X)inv(X);
% A=zeros(3,3,2);
% A(:,:,1)=eye(3);
% A(:,:,2)=magic(3);
% AInv=applyFunToEachMatrix(f,A)
%One will find that all of the matrices in A have been inverted.
%
%EXAMPLE 2:
%In this instance, we multiply a number of matrices.
% f=@(X,Y)(X*Y);
% A=zeros(3,3,2);
% A(:,:,1)=eye(3);
% A(:,:,2)=magic(3);
% B=zeros(3,3,2);
% B(:,:,1)=toeplitz([3;2;1]);
% B(:,:,2)=toeplitz([1;2;3]);
% ABProd=applyFunToEachMatrix(f,A,B)
%
%February 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

ASize=size(A);

if(length(ASize)<=2)
    numMats=1;
else
    numMats=prod(ASize(3:end));
end

if(nargin>3&&~isempty(matDims))
    X=zeros([matDims,ASize(3:end)]);
else
    X=zeros(ASize);
end

if(nargin<3||isempty(B))
    for curMat=1:numMats
        X(:,:,curMat)=f(A(:,:,curMat)); 
    end
else
    for curMat=1:numMats
        X(:,:,curMat)=f(A(:,:,curMat),B(:,:,curMat)); 
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
