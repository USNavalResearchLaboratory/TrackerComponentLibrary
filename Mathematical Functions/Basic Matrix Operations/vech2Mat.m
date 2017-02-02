function A=vech2Mat(theVec,isSymmetric,offDiagScalFact)
%%VECH2MAT Given the main diagonal and lower-triangular elements of a
%          matrix stacked column-wise, obtain the original
%          symmetric or lower-triangular matrix. This function is
%          essentially the inverse of the vech function when applied to a
%          symmetric or lower-triangular matrix with the proper
%          parameterization.
%
%INPUTS: theVect An nX1 or 1Xn vector holding the diagonal and
%                lower-triangular components of a symmetric matrix stacked
%                column-wise.
%    isSymmetric An optional boolean parameter specifying whether the
%                matrix being recreated is symmetric. If false, a lower-
%                triangular matrix is formed. The default if this parameter
%                is omitted or an empty matrix is passed is true.
%   offDiagScalFact An optional scalar factor by which the off-diagonal
%                elements are scaled. If this parameter is omitted or an
%                empty matrix is passed, then the default value of 1 is
%                used.
%
%OUTPUTS: A The matrix impled by theVec and isSymmetric. This is a dXd
%           matrix, where n=d*(d+1)/2.
%
%The opposite of this function is vech.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(offDiagScalFact))
    offDiagScalFact=1;
end

if(nargin<2||isempty(isSymmetric))
   isSymmetric=true; 
end

n=length(theVec);
theVec=theVec(:);

d=(1/2)*(-1+sqrt(1+8*n));
A=zeros(d,d);

curStart=1;
if(isSymmetric)
    for curCol=1:d
        %The diagonal element
        A(curCol,curCol)=theVec(curStart);
        span=(curStart+1):(curStart+d-curCol);
        
        %The lower-triangular part
        A((curCol+1):end,curCol)=theVec(span)*offDiagScalFact;
        %The upper triangular part.
        A(curCol,(curCol+1):end)=theVec(span)'*offDiagScalFact;
        curStart=curStart+d-curCol+1;
    end
else
    for curCol=1:d
        %The diagonal element
        A(curCol,curCol)=theVec(curStart);
        span=(curStart+1):(curStart+d-curCol);
        
        %The lower-triangular part
        A((curCol+1):end,curCol)=theVec(span)*offDiagScalFact;

        %The upper-triangular part is just zero.
        curStart=curStart+d-curCol+1;
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
