function Sk1Inv=invSymMatDelRowCol(SkInv,rowNum)
%%INVSYMMATDELROWCOL Given the inverse of a kXk real symmetric matrix Sk,
%                    one desired the inverse of the matrix obtained be
%                    removing the row/column selected by rowNum from Sk.
%                    The inverse of the reduced matrix is found using an
%                    identity that involves a single division, but not any
%                    matrix inversion.
%
%INPUTS: SkInv The kXk inverse of the real symmetric matrix Sk.
%       rowNum The index of the row/column to be deleted from Sk.
%
%OUTPUTS: Sk1Inv The (k-1)X(k-1) inverse of the real symmetric matrix that
%                is obtained by deleting row and column rowNum from Sk.
%
%This type of updating has uses in algorithms, such as in [1], where the
%formulae used are explicitly given.
%
%EXAMPLE:
%Here, we look at the relative error of this function compared to using Sk1
%and directly inverting it
% Sk=[140,   4,  37,  34,  49,  28;
%       4, 134,  16,  49,  28,  61;
%      37,  16,  74,  55,  61,  49;
%      34,  49,  55, 104,  22,  28;
%      49,  28,  61,  22,  98,  34;
%      28,  61,  49,  28,  34,  92];
% k=size(Sk,1);
% SkInv=inv(Sk);
% rowNum=4;
% Sk1Inv=invSymMatDelRowCol(SkInv,rowNum);
% sel=[1:(rowNum-1),(rowNum+1):k];
% Sk1=Sk(sel,sel);
% Sk1InvDirect=inv(Sk1);
% RelError=(Sk1InvDirect-Sk1Inv)./abs(Sk1InvDirect)
%One will see that the relative error is less than 1e-13, which is
%around what one might expect due to finite precision errors.
%
%REFERENCES:
%[1] P. M. Narendra and K. Fukunaga, "A branch and bound algorithm for
%    feature subset selection," IEEE Transactions on Computers, vol. C-26,
%    no. 9, pp. 917-922, Sep. 1977.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

k=size(SkInv,1);

sel=[1:(rowNum-1),(rowNum+1):k];
A=SkInv(sel,sel);
b=SkInv(rowNum,rowNum);
c=SkInv(sel,rowNum);

Sk1Inv=A-(1/b)*(c*c');

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
