function SkInv=invSymMatAddRowCol(Sk1Inv,y,skk)
%%INVSYMMATADDROWCOL Given the inverse of a (k-1)X(k-1) real symmetric Sk1
%               matrix, one desires the inverse of the kXk matrix that is
%               Sk1 augmented by an additional row column, where skk is the
%               new diagonal element and y is the values in the column
%               above the diagonal. The inverse of the augmented matrix is
%               found using an identity that involves a single division,
%               but not any matrix inversion. The augmented matrix
%               must be positive definite.
%
%INPUTS: SkiInv The (k-1)X(k-1) inverse of a real, symmetric matrix Sk1.
%             y The (k-1)X1 column to be added to the matrix Sk1, without
%               the diagonal element of Sk.
%           skk The diagonal element of the matrix Sk.
%
%OUTPUTS: SkInv The inverse of the kXk symmetric matrix that is Sk
%               augmented by y and skk.
%
%This type of updating has uses in algorithms, such as in [1], where the
%formulae used are explicitly given.
%
%EXAMPLE:
% Sk=[140,   4,  37,  34,  49,  28;
%       4, 134,  16,  49,  28,  61;
%      37,  16,  74,  55,  61,  49;
%      34,  49,  55, 104,  22,  28;
%      49,  28,  61,  22,  98,  34;
%      28,  61,  49,  28,  34,  92];
% k=size(Sk,1);
% Sk1=Sk(1:(k-1),1:(k-1));
% Sk1Inv=inv(Sk1);
% SkInvDirect=inv(Sk);
% y=Sk(1:(k-1),k);
% skk=Sk(end,end);
% SkInv=invSymMatAddRowCol(Sk1Inv,y,skk);
% RelError=(SkInvDirect-SkInv)./abs(SkInvDirect)
%One will see that the relative error is less than 1e-14, which is
%around what one might expect due to finite precision errors.
%
%REFERENCES:
%[1] P. M. Narendra and K. Fukunaga, "A branch and bound algorithm for
%    feature subset selection," IEEE Transactions on Computers, vol. C-26,
%    no. 9, pp. 917-922, Sep. 1977.
%
%August 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

dInv=1/(skk-y'*Sk1Inv*y);

diagTerm=-dInv*Sk1Inv*y;
prodTerm=Sk1Inv*y;

SkInv=[Sk1Inv+dInv*(prodTerm*prodTerm'), diagTerm;
	                          diagTerm', dInv];
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
