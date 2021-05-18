function theVec=vech(H,offDiagScalFact,omitDiagonal)
%%VECH The vector half operator. Given a square 2D matrix H, this returns
%      the entries that are on or below the main diagonal of H, listed
%      column-by-column, unless omitDiagonal=true, in which case the main
%      diagonal is omitted too.
%
%INPUTS: H An nXn matrix. It can be real or complex.
% offDiagScalFact An optional scalar factor by which the off-diagonal
%          elements are scaled. If this parameter is omitted or an empty
%          matrix is passed, then the default value of 1 is used.
% omitDiagonal An optional parameter indicating whether the elements on the
%          main diagonal should be omitted. The default if omitted or an
%          empty matrix is passed is false.
%
%OUTPUTS: theVec A (n*(n+1)/2)X1 matrix holding the lower-triangular part
%                of the matrix, listed column-by-column. If omitDiagonal is
%                true, then this is a (n*(n+1)/2-n)X1 vector, as the
%                elements on the main diagonal were omitted.
%
%The opposite of this function is vech2Mat.
%
%December 2015 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(omitDiagonal))
    omitDiagonal=false;
end

if(nargin<2||isempty(offDiagScalFact))
    offDiagScalFact=1;
end

n=size(H,1);

if(omitDiagonal)
    numEls=n*(n+1)/2-n;
    theVec=zeros(numEls,1);
    if(numEls==0)
        return;
    end
    curStart=1;
    for curCol=1:n
        %Off-diagonal elements
        span=(curStart):(curStart+n-curCol-1);
        theVec(span)=H((curCol+1):end,curCol)*offDiagScalFact;

        curStart=curStart+n-curCol;
    end
else
    theVec=zeros(n*(n+1)/2,1);

    curStart=1;
    for curCol=1:n
        %The diagonal element
        theVec(curStart)=H(curCol,curCol);

        %Off-diagonal elements
        span=(curStart+1):(curStart+n-curCol);
        theVec(span)=H((curCol+1):end,curCol)*offDiagScalFact;

        curStart=curStart+n-curCol+1;
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
