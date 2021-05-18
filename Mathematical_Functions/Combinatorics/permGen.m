function permVal=permGen(A,p)
%%PERMGEN Compute a generalized permanent of a rectangular matrix. For
%         numRows<numCols, the standard matrix permanent is the sum of all
%         possible products of choosing one item from each row and at most
%         one item from each column. (for a square matrix, it is the
%         determinant with all - signs changed to plus signs). In the
%         generalized permanet here, p(1) things are chosen from the first
%         column, p(2), from the second column, etc. (sum(p)<=number of
%         rows). No row is chosen twice. All chosen values are multiplied
%         together and then summed over all possible choices. The
%         generalized matrix permanent arises when computing
%         multihomogenerous Bézout bounds.
%
%INPUTS: A An nXm matrix with m>=n.
%        p A mX1 vector of positive integer values such that each p(i)>=1
%          and sum(p)=n. This is a list of the number of things to take
%          from each row for the generalized permanent.
%
%OUTPUTS: permVal The value of the generalized matrix permanent.
%
%The generalized permanent of a rectangular matrix arises in the
%computation of the multihomogenerous Bézout bound in Chapter 3 of [1].
%Unlike the standard permanent, there does not appear to be a particularly
%efficient algorithm in the literature. Thus, this function is just a
%brute-force implementation, generating every term of the sum.
%
%EXAMPLE:
%As in the multihomogenerous Bézout bound in Section 2.2. of [1], we have
% A=[2,0,0;
%    0,2,0;
%    0,0,2;
%    0,2,0;
%    1,1,0;
%    1,1,0;
%    1,1,1;
%    1,0,0];
% p=[2;4;2];
% permVal=permGen(A,p)
%One will find permVal=16.
%
%REFERENCES:
%[1] J. Verschelde, "Homotopy continuation methods for solving polynomial
%    systems," Ph.D. dissertation, Katholieke Universiteit Leuven, Leuven,
%    Belgium, May 1996.
%
%March 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

numRows=size(A,1);
numCols=size(A,2);

if(sum(p)~=numRows||any(p<=0)||~isreal(p))
   error('The elements in p must be >=1 and sum to the number of rows of the matrix.') 
end

p=p(:);

%If the number of rows equals the number of columns, then it is the same as
%the standard matrix permanent algorithm and the perm function is faster.
if(numRows==numCols)
    permVal=perm(A);
    return;
end

permVal=0;

curColumn=1;
descending=true;
prodVals=zeros(numCols,1);

%To hold the indices of rows that can be assigned in the current column.
freeRows=zeros(numRows,numCols);
freeRows(:,1)=1:numRows;

%The number of rows that are left for each column.
freeRowNum=[numRows;numRows-cumsum(p(1:(end-1)))];

%Holds the current combinations of assigned rows.
assignedCombo=zeros(max(p),numCols);

while(curColumn>0)
    if(descending)
        %Just entering this column.
        if(curColumn==numCols)
            %If at the final column, loop though all the possible combos
            %and add the final product values to permVal.
            prevProd=prod(prodVals(1:(curColumn-1)));
            
            %The first combination.
            curCombo=0:(p(curColumn)-1);
            %Go through all of the combinations.
            while(~isempty(curCombo))
                assignedIdx=curCombo+1;
                assignedRows=freeRows(assignedIdx,curColumn);
                prodVal=prod(A(assignedRows,curColumn));
                
                permVal=permVal+prevProd*prodVal;
                curCombo=getNextCombo(curCombo,freeRowNum(curColumn));
            end

            descending=false;
            curColumn=curColumn-1;
            continue;
        else
            %The first combination of assigned rows for this column.
            %assignedCombo is actually the indices of the current column of
            %freeRowNum that are assigned.
            assignedCombo(1:p(curColumn),curColumn)=0:(p(curColumn)-1);
        end  
    else%Going up
        %Advance to the next combination.
        newCombo=getNextCombo(assignedCombo(1:p(curColumn),curColumn),freeRowNum(curColumn));
        if(isempty(newCombo))%Keep going up if last combination reached.
            curColumn=curColumn-1;
            continue;
        end
        %Otherwise, assign the combination and stop descending.
        assignedCombo(1:p(curColumn),curColumn)=newCombo;
        descending=true;
    end
    
    assignedIdx=assignedCombo(1:p(curColumn),curColumn)+1;
    assignedRows=freeRows(assignedIdx,curColumn);
    prodVals(curColumn)=prod(A(assignedRows,curColumn));
    
    %The final result would hvae just been zero, so do not bother with
    %increasing the columns anymore. Just act as if coming up.
    if(prodVals(curColumn)==0)
        descending=false;
        continue;
    end
    
    %Fill in the row numbers in freeRowNum for the next column.
    sel=true(freeRowNum(curColumn),1);
    sel(assignedIdx)=false;
    freeRows(1:freeRowNum(curColumn+1),curColumn+1)=freeRows(sel,curColumn);
    curColumn=curColumn+1;
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
