function [tuples,gain]=assign2DFreeCol(C,maximize)
%%ASSIGN2DFREECOL Solve the two-dimensional assignment problem where the
%          binary constraints on the columns have been removed. That is,
%          minimize (or maximize)
%          \sum_{i=1}^{numRow}\sum_{j=1}^{numCol}C_{i,j}*x_{i,j}
%          subject to 
%          \sum_{j=1}^{numCol}x_{i,j} =1 for all i
%          x_{i,j}=0 or 1
%          This problem is just a naïve nearest neighbor assignment; one
%          simply assigns each row to the column with the smallest (or
%          largest) value.
%
%INPUTS: C A numRowXnumCol cost matrix that does not contain any NaNs.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          or an empty matrix is passed is false.
%
%OUTPUTS: tuples A 2XnumRow set of assignment values. This is ordered
%                [rowIndex;columnIndex].
%          gain This is the value of the cost. This is the sum of the
%               values in C corresponding to the tuples.
%
%November 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<2||isempty(maximize))
    maximize=false;
end

numRow=size(C,1);
tuples=zeros(2,numRow);
gain=0;

if(maximize)
    for curRow=1:numRow
        [val,maxIdx]=max(C(curRow,:));

        tuples(:,curRow)=[curRow;maxIdx];
        gain=gain+val;
    end
else
    for curRow=1:numRow
        [val,minIdx]=min(C(curRow,:));

        tuples(:,curRow)=[curRow;minIdx];
        gain=gain+val;
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
