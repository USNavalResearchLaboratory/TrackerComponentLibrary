function A=triLower(A)
%%TRILOWER Zero the upper-triangular potion of the matrix and multiply the
%          diagonal terms by 0.5
%
%INPUTS: A A square matrix.
%
%OUTPUTS:  A The matrix A with the upper-triangular portion zeroed and the
%            diagonal elements divided by 2.
%
%This function plays a role in turning the derivative of the covariance
%matrix into the derivative of the lower-triangular square root of the
%covariance matrix. It is described in Equation (80) in [1].
%
%REFERENCES
%[1] D. F. Crouse, "Basic tracking using nonlinear 3D monostatic and
%    bistatic measurements," IEEE Aerospace and Electronic Systems 
%    Magazine, vol. 29, no. 8, Part II, pp. 4-53, Aug. 2014.
%
%October 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numRows=size(A,1);
numCol=size(A,2);

for curRow=1:numRows
    A(curRow,curRow)=0.5*A(curRow,curRow);
    A(curRow,(curRow+1):numCol )=0;
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
