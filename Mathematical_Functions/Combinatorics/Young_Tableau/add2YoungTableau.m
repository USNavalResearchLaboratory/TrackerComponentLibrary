function [P,s,t]=add2YoungTableau(x,P)
%%ADD2YOUNGTABLEAU Given a Young tableau of integers, add an integer value
%           to the tableau that is not already in the tableau. A Young
%           tableau is an arrangement of integers in m rows such that the
%           number of items in each row is ordered n1>=n2>=n3>=...>=nm,
%           the elements in each row are left-justified and the elements in
%           each row are given in increasing order. Young tableaus are
%           related to involutions (permutations that are their own
%           inverses). Each Young tableau can be transformed into an
%           equivalent involution.
%
%INPUTS: x The integer >=1 that is to be added to the Young tableau. THis
%          should not already be in P.
%        P An MXN matrix holding the Young tableau. Elements that are not
%          part of the tableau should be set to Inf. The matrix can have
%          more rows and columns than are necessary for the tablau, with
%          all extra elements set to Inf. Omitting P or passing an empty
%          matrix just starts a new tableau.
%
%OUTPUTS: P The tableau P with x added. The size of P compared to the input
%           may have increased and there may be an extra row or column with
%           Inf elements. If P on the input had at least one extra row and
%           at least one extra column, then its size is guaranteed not to
%           increase in either dimension.
%       s,t The indices where a value that was previously not in the
%           tableau (e.g might be Inf or outside of the dimensions of P)
%           was added.
%
%This function implements Algorithm I in Section 5.1.4 of [1].
%
%EXAMPLE:
%Here, we implement the example prior to Algorithm I in Section 5.1.4 of
%[1].
% P=[1,   3,  5,  9,  12, 16;
%    2,   6,  10, 15, Inf,Inf;
%    4,   13, 14, Inf,Inf,Inf;
%    11,  Inf,Inf,Inf,Inf,Inf;
%    17,  Inf,Inf,Inf,Inf,Inf];
% [P,s,t]=add2YoungTableau(8,P)
%One gets
% P =[ 1,   3,   5,   8,  12,  16, Inf;
%      2,   6,   9,  15, Inf, Inf, Inf;
%      4,  10,  14, Inf, Inf, Inf, Inf;
%     11,  13, Inf, Inf, Inf, Inf, Inf;
%     17, Inf, Inf, Inf, Inf, Inf, Inf;
%    Inf, Inf, Inf, Inf, Inf, Inf, Inf];
% and (s,t)=(4,2)
%
%REFERENCES:
%[1] D. Knuth, The Art of Computer Programming: Sorting and Searching, 2nd
%    ed. Boston: Addison-Wesley, 2018, vol. 3.
%
%November 2018 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%If starting a new tableau.
if(nargin<2||isempty(P))
    P=x;
    s=1;
    t=1;
    return;
elseif(all(~isfinite(P(:))))
    P(1,1)=x;
    s=1;
    t=1;
    return;
end

numRow=size(P,1);
numCol=size(P,2);

augmentCol=isfinite(P(1,numCol));
augmentRow=isfinite(P(numRow,1));

%If P must be increased.
if(augmentRow||augmentCol)
    PNew=Inf*ones(numRow+augmentRow,numCol+augmentCol);
    PNew(1:numRow,1:numCol)=P;
    P=PNew;
    numRow=numRow+augmentRow;
end

%Allocate space.
xi=Inf*ones(numRow+1,1);

%Step I1
i=1;
xi(1)=x;
j=find(P(1,:)==Inf,1);

while(1)
    %Step I2, find xi(i+1)
    while(j>1&&xi(i)<P(i,j-1))
        j=j-1;
    end
    xi(i+1)=P(i,j);

    %Step I3
    P(i,j)=xi(i);

    %Step I4
    if(isfinite(xi(i+1)))
        i=i+1; 
        continue;
    else
        %Step I5
        s=i;
        t=j;
        return
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
