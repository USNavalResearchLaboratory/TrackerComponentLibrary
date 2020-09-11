function CNew=invAssign2D(C,col4rowDes,maximize,epsVal)
%%INVASSIGN2D Solve the inverse 2D assignment problem. That is, given a
%             cost matrix and a feasible but suboptimal solution, make the
%             minimal changes to the cost matrix possible such that the
%             suboptimal solution becomes optimal. We are minimizing
%             sum(abs(C(:)-CNew(:))), when only considering the finite
%             elements in C. This function is made for square cost
%             matrices.
%
%INPUTS:C An NXN cost matrix that does not contain any NaNs. Forbidden
%          assignments can be given costs of +Inf for minimization and -Inf
%          for maximization.
% col4rowDes An NX1 or 1XN vector giving the column assigned to each row in
%          the desired optimal solution. This vector must be a feasible
%          solution.
% maximize If true, the minimization problem is transformed into a
%          maximization problem. The default if this parameter is omitted
%          is false.
%   epsVal In general, this function will only make the desired solution
%          have the same cost as the current optimal solution.
%          Additionally, given finite precision limitations, it is also
%          possible that the desired optimal solution will be suboptimal by
%          a tiny epsilon value. This value is what is subtracted from the
%          assignment costs to help guarantee that the desired solution is
%          optimal (and not off by some epsilon). If omitted or an empty
%          matrix is passed, this is 2^8*eps(val) where val is the maximum
%          finite absolute value of the elements in C. Generally, that
%          should suffice. However, one might want to make this value zero
%          when dealing with all-integer costs where one knows that finite-
%          precision effects will not be an issue.
%
%OUTPUTS: CNew The modified NXN cost matrix where col4rowDes is an optimal
%              solution.
%
%The algorithm is taken from Section 4 of [1].
%
%EXAMPLE 1:
%This is Example 4.4 of [1].
% C=[5,7,3,Inf;
%    4,7,4,2;
%    Inf,3,8,3;
%    8,7,4,6];
% col4rowDes=[1;2;3;4];
% %First, we see that col4RowDes is not currently the optimal solution.
% col4row=assign2D(C)
% %One will see that col4row=[1;4;2;3]. Now, we modify the cost matrix.
% CNew=invAssign2D(C,col4rowDes);
% %The new optimal assignment is col4rowDes:
% col4rowNew=assign2D(CNew);
% all(col4rowNew==col4rowDes)
%The result is true (1).
%
%EXAMPLE 2:
%This example shows that while the inverse function DOES indeed make the
%desired solution optimal, it does NOT guarantee that it will be unique if
%epsVal=0. This is an example of where running assign2D after invAssign2D,
%one gets an optimal assignment DIFFERENT than the desired solution.
%However, after comparing, one will see that both solutions have the same
%gain and are thus both optimal. Note, however, that with non-integer
%values for C and epsVal=0, it is possible that the deisred solution could
%be marginally suboptimal (by a value within a few magnitudes of eps() of
%the largest finite magnitude value in C).
% C=[-36     4    71   -76   -68    18;
%     70   174    31   -81   -66    55;
%    141    15    41    51    86    68;
%   -160  -123   -57    -1    11   117;
%    102  -219    14  -115    39    47;
%    145   -33  -163     0    88   141];
% N=length(C);
% col4rowDes=[4;2;5;3;6;1];
% 
% %Modify the cost matrix.
% CNew=invAssign2D(C,col4rowDes,false,0);
% %The new optimal assignment is col4rowDes:
% col4rowAlt=assign2D(CNew);
% 
% gainValDes=0;
% gainValAlt=0;
% for k=1:N
%     gainValDes=gainValDes+CNew(k,col4rowDes(k));
%     gainValAlt=gainValAlt+CNew(k,col4rowAlt(k));
% end
% gainValDes==gainValAlt
% all(col4rowDes==col4rowAlt)
%One will see that the gains are equal, but the assignments are different.
%
%REFERENCES:
%[1] J. Zhang and Z. Liu, "Calculating some inverse linear programming
%    problems," Journal of Computational and Applied Mathematics, vol.72,
%    no. 2, pp. 261-273, 13 Aug. 1996.
%
%August 2019 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3||isempty(maximize))
    maximize=false;
end

numRow=size(C,1);
numCol=size(C,2);

if(numRow~=numCol)
  error('This algorithm requires that numRow==numCol')
end

if(nargin<4||isempty(epsVal))
    maxVal=max(abs(C(isfinite(C))));
    epsVal=2^8*eps(maxVal);
end

%The cost matrix must have all non-negative elements for the assignment
%algorithm to work. This forces all of the elements to be positive. The
%delta is added back in when computing the gain in the end. Positiveity is
%a requirement of the algorithm used in the assign2D function. If we do not
%make this adjustment prior to calling assign2D, then assign2D will make it
%and the dual variables returned will not be suitable for adjusting the
%matrix.
if(maximize==true)
    CDelta=max(C(:));

    %If C is all negative, do not shift.
    if(CDelta<0)
        CDelta=0;
    end

    C=-C+CDelta;
else
    CDelta=min(C(:));

    %If C is all positive, do not shift.
    if(CDelta>0)
        CDelta=0;
    end

    C=C-CDelta;
end

%Step 1 after Theorem 4.3 in [1].
%Rearrange the columns so that the desired optimal solution will be on the
%diagonals.
C=C(:,col4rowDes);

%Solve the 2D assignment problem on the cost matrix with the rearranged
%entries. The definitions of u and v are swithed in [1] compared to
%assign2D.
[~,~,~,v,u]=assign2D(C,maximize);

CNew=C;
for curRow=1:numRow
    CNew(curRow,curRow)=u(curRow)+v(curRow)-epsVal;
end

%Determine the reverse permutation to put the columns back in the orignal
%order.
invPerm=inversePermutation(col4rowDes);

%Undo the permutation at the beginning.
CNew=CNew(:,invPerm);

%Undo the original transformation of the matrix.
if(maximize==true)
    CNew=-(CNew-CDelta);
else
    CNew=CNew+CDelta;
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
