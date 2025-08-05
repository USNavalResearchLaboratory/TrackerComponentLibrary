function boolVals=pointIsInEllipsoid(xPts,xCen,A,matType,gamma)
%%POINTISINELLIPSOID Determine whether or not a set of points are within an
%       ellipsoid. A point x is in the ellipsoid if
%       (x-xCen)'*inv(A)*(x-xCen)<=gamma if Type=0 or
%       (x-xCen)'*inv(A*A')*(x-xCen)<=gamma if Type=1 or
%       (x-xCen)'*A*(x-xCen)<=gamma if Type=2.
%
%INPUTS: xPts A numDimXN set of points to test.
%     xCen The numDimX1 center of the ellipsoid.
%        A The numDimXnumDim matrix defining the shape fo the ellipsoid.
%          The form of A is determined by MatType.
%  matType This optional parameter specified the type of matrix that M is.
%          Possible values are
%          0 (The default if omitted or an empty matrix is passed) A is the
%            matrix A in the form (x-xCen)'*inv(A)*(x-xCen)<=gamma.
%          1 A is the matrix A in the form
%            (x-xCen)'*inv(A*A')*(x-xCen)<=gamma.
%          2 A is the matrix A in the form (x-xCen)'*A*(x-xCen)<=gamma.
%    gamma A positive real threshold. The default if omitted or an empty
%          matrix is passed is zero.
%
%OUTPUTS: boolVals A 1XN set of boolean values indicating whether or not
%                  the point is in the ellipsoid.
%
%The function calls invSymQuadForm for matType 0 and 1.
%
%June 2023 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<5||isempty(gamma))
    gamma=1;
end

if(nargin<4||isempty(matType))
    matType=0;
end

if(matType==2)
    N=size(xPts,2);
    vals=zeros(1,N);
    for k=1:N
        diffVal=xPts(:,k)-xCen;
        vals(k)=diffVal'*A*diffVal;
    end
else
    diffVals=bsxfun(@minus,xPts,xCen);
    vals=invSymQuadForm(diffVals,A,matType);
end

boolVals=vals<=gamma;

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
