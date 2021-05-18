function C=circulantMatrix(v,reverseDir)
%%CIRCULANTMATRIX Create a circulant matrix where the first row is given by
%            v. A circulatn matrix is a matrix where each row going down
%            from the top rotates the elements of the previous row by 1.
%
%INPUTS: v A 1XN or NX1 vector containing the elemtns that go into the
%          first row.
% reverseDir If true, this boolean parameter indicates that each subsequent
%          row should be left-shifted from the previous one. If this
%          parameter is omitted or an empty matrix is passed, then the
%          default is false.
%
%OUTPUTS: C The NXN circulant matrix whose first row is v.
%
%EXAMPLE:
% CF=circulantMatrix(1:5,false)
% CR=circulantMatrix(1:5,true)
%The results are:
%CF=[1,2,3,4,5;
%    5,1,2,3,4;
%    4,5,1,2,3;
%    3,4,5,1,2;
%    2,3,4,5,1];
%CR=[1,2,3,4,5;
%    2,3,4,5,1;
%    3,4,5,1,2;
%    4,5,1,2,3;
%    5,1,2,3,4];
%
%October 2017 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

    if(nargin<2||isempty(reverseDir))
        reverseDir=false;
    end

    N=length(v);
    v=v(:);
    x=[v(N:-1:1,1);flipud(v(2:end))];
    if(reverseDir)
        idx=bsxfun(@plus,mod(0:-1:-((N-1)),N)',(N:-1:1));
    else
        idx=bsxfun(@plus,(0:(N-1))',(N:-1:1));
    end
    
    C=x(idx);
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
