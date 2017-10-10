function z=conv2Sparse(x,y,shape)
%%CONV2SPARSE Compute the 2D convolution of two two-dimensional sparse
%             matrices. This is the same as the conv2 function, except x
%             and y are sparse. A 2D convolution is given (with indexation
%             from zero, not 1) as
%     z(i,j)=sum_{m=0}^{numRowsX-1}sum_{n=0}^{numRowsY-1} x(m,n)*y(i-m,j-n)
%             for 0<=i<numRowsX+numRowsY-1 and 0<=j<numRowsX+numRowsY-1.
%             This is the "full" version of the convolution. The 'valid'
%             and 'same' variants described below discard certain elements
%             of the output. 2D convolutions are often used in image
%             processing. Sparse variants tend to arise more in centroiding
%             algorithms.
%
%INPUTS: x,y A numRowsXXnumColsX and a numRowsYXnumColsY pair of sparse
%            matrices. These can be complex.
%      shape An optional parameter that determines the shape of the output
%            matrix x. Possible values are
%            'full' (The default if omitted or an empty matrix is passed).
%                    Return the full convolution.
%             'same' Return only the central portion of the convolution.
%                    The output is the same size as x.
%            'valid' Only return the part of the convolution that has no
%                    zero padding. When moving matrix y over x (as the sum
%                    in the convolution can be viewed for different i and
%                    j), this is essentially getting rid of times where any
%                    elements in y go beyond the end of x.
%
%OUTPUTS: z A sparse matrix that holds the 2D convolution of x and y. The
%           output should be the same as conv2(full(x),full(y),shape).
%
%January 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(shape))
    shape='full'; 
end

[rowX,colX,valX]=find(x);
[rowY,colY,valY]=find(y);
%Note that the outputs of ndgrid with two inputs are transposes of the
%outputs of meshgrid.
[rowGrid1,rowGrid2]=ndgrid(rowX,rowY);
[colGrid1,colGrid2]=ndgrid(colX,colY);
allProdVals=valX(:)*valY(:).';

[numRowX,numColX]=size(x);
[numRowY,numColY]=size(y);
switch(shape)
    case 'full'
        %numRowX+numRowY-1 rows and numColX+numColY-1 columns in the full
        %convolution.
        numElsX=numRowX+numRowY-1;
        numElsY=numColX+numColY-1;
        
        rowIdx=rowGrid1(:)+rowGrid2(:)-1;
        colIdx=colGrid1(:)+colGrid2(:)-1;
    case 'valid'
        numElsX=max(0,numRowX-max(0,numRowY-1));
        numElsY=max(0,numColX-max(0,numColY-1));

        rowIdx=rowGrid1(:)+rowGrid2(:)-numRowY;
        colIdx=colGrid1(:)+colGrid2(:)-numColY;
        
        sel=(rowIdx>=1)&(rowIdx<=numElsX)&(colIdx>=1)&(colIdx<=numElsY);
        rowIdx=rowIdx(sel);
        colIdx=colIdx(sel);
        allProdVals=allProdVals(sel);
    case 'same'
        numElsX=numRowX;
        numElsY=numColX;
        
        rowIdx=rowGrid1(:)+rowGrid2(:)-ceil((numRowY+1)/2);
        colIdx=colGrid1(:)+colGrid2(:)-ceil((numColY+1)/2);
        
        sel=(rowIdx>=1)&(rowIdx<=numRowX)&(colIdx>=1)&(colIdx<=numColX);
        rowIdx=rowIdx(sel);
        colIdx=colIdx(sel);
        allProdVals=allProdVals(sel);
    otherwise
        error('Unknown shape specified')
end

z=sparse(rowIdx,colIdx,allProdVals(:),numElsX,numElsY);

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
