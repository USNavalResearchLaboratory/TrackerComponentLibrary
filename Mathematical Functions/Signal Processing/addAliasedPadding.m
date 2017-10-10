function XPadded=addAliasedPadding(X,numPad,dim)
%%ADDALIASEDPADDING Increase the size of the one-dimensional array or two-
%           dimensional matrix X by padding. For a matrix, pad numPad(1)
%           elements before the rows start and after they end and numPad(2)
%           elements before the columns start and after they end.
%           Alternatively only rows or only columns can be padded. This
%           padding is filled with the values of X aliased back on itself.
%           This type of padding can be useful when performing cell
%           averaging CFAR or when centroiding.
%
%INPUTS: X The numRowsXnumColsXnumX set of numX matrices to which aliased
%          padding should be added.
%   numPad A 2X1 vector where numPad(1) is the amount to pad the rows on
%          each side and numPad(2) is the amount to pad the columns on each
%          side. If only one dimension is to be padded, then the dimension
%          padded in given by dim and this value is a scalar. It is
%          required that numPad(1)<=numRows and numPad(2)<=numCols when
%          padding only rows or only columns. If a scalar is passed and
%          both dimensions are to be padded, then it is assumed that the
%          same amount of padding is to be used on both dimensions.
%      dim This specifies the dimension to pad when numPad is scalar.
%          dim=1 pads the rows; dim=2 pads the columns and omitting dim,
%          passing an empty matrix, or setting dim=0 pads the rows and
%          columns by the same amount. The default if this parameter is
%          omitted or an empty matrix is passed is dim=0, unless X is a
%          vector in which case the non-singletom dimension is the one that
%          is padded.
%
%OUTPUTS: XPadded The (numRows+2*numPad(1))X(numRows+2*numPad(2))XnumX set
%                 of numX matrix where the central portion is X and all of
%                 the sides have values that have been aliased from X.
%
%EXAMPLE:
% X=magic(4);
% numPad=2;
% XPadded=addAliasedPadding(X,numPad)
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numRows=size(X,1);
numCols=size(X,2);
numX=size(X,3);

if(nargin<3||isempty(dim))
    if(numRows>1&&numCols>1||numRows==1&&numCols==1)
        dim=0;
    elseif(numRows>1&&numCols==1)
        dim=1;
    else%(numRows==1&&numCols>1)
        dim=2;
    end
end

if(length(numPad)==1)
    numPad=[numPad;numPad];
end

switch(dim)
    case 0%Pad both
        XPadded=zeros(numRows+2*numPad(1),numCols+2*numPad(2),numX);
        
        %Indices of elements that will be selected in XPadded.
        innerRows=(numPad(1)+1):(numPad(1)+numRows);
        innerCols=(numPad(2)+1):(numPad(2)+numCols);
        outerRowsU=1:numPad(1);
        outerRowsD=(numPad(1)+numRows+1):(numRows+2*numPad(1));
        outerColsL=1:numPad(2);
        outerColsR=(numPad(2)+numCols+1):(numCols+2*numPad(2));

        %Indices of elements that will be selected in X.
        rowPadU=(numRows-numPad(1)+1):numRows;
        rowPadD=1:numPad(1);
        colPadL=(numCols-numPad(2)+1):numCols;
        colPadR=1:numPad(2);

        %Fill in each section of XPadded.
        XPadded(innerRows,innerCols,:)=X;
        XPadded(innerRows,outerColsL,:)=X(:,colPadL,:);
        XPadded(innerRows,outerColsR,:)=X(:,colPadR,:);
        XPadded(outerRowsU,innerCols,:)=X(rowPadU,:,:);
        XPadded(outerRowsD,innerCols,:)=X(rowPadD,:,:);
        XPadded(outerRowsU,outerColsL,:)=X(rowPadU,colPadL,:);
        XPadded(outerRowsD,outerColsL,:)=X(rowPadD,colPadL,:);
        XPadded(outerRowsD,outerColsR,:)=X(rowPadD,colPadR,:);
        XPadded(outerRowsU,outerColsR,:)=X(rowPadU,colPadR,:);
    case 1%Pad the rows
        XPadded=zeros(numRows+2*numPad(1),numCols,numX);
        
        innerRows=(numPad(1)+1):(numPad(1)+numRows);
        outerRowsU=1:numPad(1);
        outerRowsD=(numPad(1)+numRows+1):(numRows+2*numPad(1));
        
        rowPadU=(numRows-numPad(1)+1):numRows;
        rowPadD=1:numPad(1);
        
        XPadded(innerRows,:,:)=X;
        XPadded(outerRowsU,:,:)=X(rowPadU,:,:);
        XPadded(outerRowsD,:,:)=X(rowPadD,:,:);
    case 2%Pad the columns
        XPadded=zeros(numRows,numCols+2*numPad(2),numX);
        
        innerCols=(numPad(2)+1):(numPad(2)+numCols);
        outerColsL=1:numPad(2);
        outerColsR=(numPad(2)+numCols+1):(numCols+2*numPad(2));
        
        colPadL=(numCols-numPad(2)+1):numCols;
        colPadR=1:numPad(2);
        
        XPadded(:,innerCols,:)=X;
        XPadded(:,outerColsL,:)=X(:,colPadL,:);
        XPadded(:,outerColsR,:)=X(:,colPadR,:);
    otherwise
        error('Invalid value of dim specified')
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
